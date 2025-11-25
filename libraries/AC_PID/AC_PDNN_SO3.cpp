/// @file	AC_PDNN_SO3.cpp
/// @brief	Generic PDNN algorithm

#include <AP_Math/AP_Math.h>
#include "AC_PDNN_SO3.h"


const AP_Param::GroupInfo AC_PDNN_SO3::var_info[] = {
    // @Param: kR
    // @DisplayName: PID Proportional Gain
    // @Description: P Gain which produces an output value that is proportional to the current error value
    AP_GROUPINFO_FLAGS_DEFAULT_POINTER("kR_xy",    0, AC_PDNN_SO3, _kR, default_kR),
    // @Param: KOmega
    // @DisplayName: PID Derivative Gain
    // @Description: D Gain which produces an output that is proportional to the rate of change of the error
    AP_GROUPINFO_FLAGS_DEFAULT_POINTER("kOmega_xy",    1, AC_PDNN_SO3, _kOmega, default_kOmega),

    AP_GROUPINFO_FLAGS_DEFAULT_POINTER("kR_z",    2, AC_PDNN_SO3, _kR_z, default_kR_z),
    
    AP_GROUPINFO_FLAGS_DEFAULT_POINTER("kOmega_z",    3, AC_PDNN_SO3, _kOmega_z, default_kOmega_z),

    AP_GROUPEND
};

// Constructor 构造函数
AC_PDNN_SO3::AC_PDNN_SO3(float initial_kR, float initial_kOmega, float initial_kR_z, float initial_kOmega_z) :
    default_kR(initial_kR),
    default_kOmega(initial_kOmega),
    default_kR_z(initial_kR_z),
    default_kOmega_z(initial_kOmega_z)
{
    // load parameter values from eeprom
    AP_Param::setup_object_defaults(this, var_info); //读取eeprom存储参数值，也可以不用eeprom，选择在代码中直接定义硬编码参数值，坏处是调试后每次都会重置，不会保存。
    
    // reset input filter to first value received 重置控制器
    _reset = true; //每次调用重置为true
}

//  update_all - set target and measured inputs to PDNN controller and calculate outputs
//  target and error are filtered
//  the derivative is then calculated and filtered
//  the integral is then updated if it does not increase in the direction of the limit vector
Vector3f AC_PDNN_SO3::update_all(const Matrix3f &R_c, const Matrix3f &R, const Vector3f &Omega, float dt, bool Rc_active)
{
    // don't process inf or NaN //检查输入的有效性，避免处理空值NaN与无穷大inf的值
    if (R_c.is_nan() || R.is_nan()) {
        return Vector3f{}; //返回一个vector3f避免报错
    }
    _R = R;
    _Omega = Omega;
    _Rc_active = Rc_active; //检查位置控制环是否被调用
    //更新_Omega_hat斜对称矩阵
    _Omega_hat.a.x = 0.0f;_Omega_hat.a.y = -_Omega.z;_Omega_hat.a.z = _Omega.y;
    _Omega_hat.b.x = _Omega.z;_Omega_hat.b.y = 0.0f;_Omega_hat.b.z = -_Omega.x;
    _Omega_hat.c.x = -_Omega.y;_Omega_hat.c.y = _Omega.x;_Omega_hat.c.z = 0.0f;
 

    // reset input filter to value received //无人机重启pdnn姿态控制时的初始化
    if (_reset) { //初始化逻辑
        _reset = false;
        
        _R_c = R_c; //更新当前循环的_R_c
        
        //旋转矩阵误差_e_R初始化
        _e_R_hat = (_R_c.transposed() * _R - _R.transposed() * _R_c) * 0.5f; //计算旋转矩阵误差的斜对称矩阵, 注意要把0.5f放在Matrix3f后面，因为函数重载的格式要求
        _e_R.x = -_e_R_hat.b.z; //斜对称矩阵.V逆运算，对_e_R进行赋值得到旋转矩阵误差_e_R，注意是一个Vector3f
        _e_R.y = _e_R_hat.a.z;
        _e_R.z = _e_R_hat.b.x;

        //Psi_R姿态误差标量函数初始化
        _Psi_R = (1.0f-(_R_c.transposed() * _R).a.x + 1.0f - (_R_c.transposed() * _R).b.y + 1.0f - (_R_c.transposed() * _R).c.z) * 0.5f;

        //防止微分爆炸，初始化微分项
        _dot_R_c.zero();
        _dot_Omega_c.zero();
        _Omega_c.zero();

        //角速度误差_e_Omega初始化
        //！！！先尝试初始化为Omega，因为前面初始化了微分项！！这里可能需要修改
        _e_Omega = _Omega;

        
        //初始化归零控制器输出
        _pdnn_output.x = 0;
        _pdnn_output.y = 0;
        _pdnn_output.z = 0;


    } else { //更新循环
        Matrix3f _R_c_last{_R_c}; //将上一个循环的_R_c存储到一个临时变量 _R_c_last 中，用于后续的微分项计算。这里用到拷贝函数，等价于Matrix3f error_last = _error;
        
        //更新当前循环的_R_c
        _R_c = R_c; 

        //更新当前循环旋转矩阵误差_e_R
        _e_R_hat = (_R_c.transposed() * _R - _R.transposed() * _R_c) * 0.5f; //计算旋转矩阵误差的斜对称矩阵，注意要把0.5f放在Matrix3f后面，因为函数重载的格式要求
        _e_R.x = -_e_R_hat.b.z; //斜对称矩阵.V逆运算，对_e_R进行赋值得到旋转矩阵误差_e_R，注意是一个Vector3f
        _e_R.y = _e_R_hat.a.z;
        _e_R.z = _e_R_hat.b.x;
       
        //更新Psi_R姿态误差标量函数
        _Psi_R = (1.0f-(_R_c.transposed() * _R).a.x + 1.0f - (_R_c.transposed() * _R).b.y + 1.0f - (_R_c.transposed() * _R).c.z) * 0.5f;

        //计算_R_c微分项，这里暂时不考虑滤波
        if (is_positive(dt)) { //检查时间步长是否有效
            _dot_R_c = (_R_c - _R_c_last) / dt;  //理论上应该可以实现逐元素求导（考虑进行正交化或者转化为四元数后归一化！！！）
        }
        
        Vector3f _Omega_c_last{_Omega_c}; //将上一个循环的_Omega_c存储到一个临时变量 _Omegac_c_last 中

        //更新当前循环的期望角速度_Omega_c
        _Omega_c_hat = _R_c.transposed() * _dot_R_c;
        _Omega_c.x = -_Omega_c_hat.b.z; //斜对称矩阵.V逆运算
        _Omega_c.y = _Omega_c_hat.a.z;
        _Omega_c.z = _Omega_c_hat.b.x;
    
        //更新当前循环角速度误差_e_Omega ！！考虑加上滤波！！！！（暂时没加）
        _e_Omega = _Omega - _R.transposed() * _R_c * _Omega_c;
        //_e_Omega = _Omega; //暂时第二项设置为0 

        //计算_Omega_c微分项，这里进行了滤波来消除微分爆炸（重要）
        if (is_positive(dt)) { //检查时间步长是否有效
            const Vector3f dot_Omega_c{(_Omega_c - _Omega_c_last) / dt}; //_error - error_last：计算当前误差与上一时刻误差之间的差，表示误差的变化量。计算误差变化量除以时间步长 dt，得到误差变化的速率，即微分项。
            _dot_Omega_c += (dot_Omega_c - _dot_Omega_c) * get_filt_D_alpha(dt);
            //_dot_Omega_c = (_Omega_c - _Omega_c_last) / dt;  //理论上应该可以实现逐元素求导，这里是不滤波的代码
        }

        //update I term 更新积分项
        //void AC_PDNN_3D::update_i(float dt, float _ki, float _c1, float _kimax, bool limit)
        update_i(dt, 1.0f, 15.0f, 30.0f, true); //尽量小，姿态控制要求实时性

       

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    }
   
   
    _pdnn_output_R.x = -_e_R.x * _kR; //计算P项输出
    _pdnn_output_R.y = -_e_R.y * _kR;
    _pdnn_output_R.z = -_e_R.z * _kR_z;

    _pdnn_output_Omega.x = -_e_Omega.x * _kOmega; //计算D项输出
    _pdnn_output_Omega.y = -_e_Omega.y * _kOmega;
    _pdnn_output_Omega.z = -_e_Omega.z * _kOmega_z;
    
    //计算几何控制项，这里惯性张量
    Matrix3f J;
    J.a.x=0.01f;J.a.y=0.0f;     J.a.z=0.0f;
    J.b.x=0.0f;      J.b.y=0.02f;J.b.z=0.0f;
    J.c.x=0.0f;      J.c.y=0.0f;     J.c.z=0.02f;
    Matrix3f J_inv;
    J_inv = J;
    J_inv.a.x = 1/J.a.x;J_inv.b.y = 1/J.b.y;J_inv.c.z = 1/J.c.z;
    //额外增广项，以实现UPAS
    Vector3f Aug;
    Aug = J_inv * _Omega_hat * J * _Omega;

    _geomrtry_output = _Omega_hat * _R.transposed() * _R_c * _Omega_c  - _R.transposed() * _R_c * _dot_Omega_c;

    //(void)_geomrtry_output;
   if (_Rc_active) {
    //计算总输出，每个方向上乘以惯性张量
    //_pdnn_output.x = _J_x * (-_e_R.x * 40.0f - _e_Omega.x * 80.0f - 0.0f * _integrator.x - _geomrtry_output.x - 1.0f *_phi_x + 0.0f*Aug.x); 
    //_pdnn_output.y = _J_y * (-_e_R.y * 40.0f - _e_Omega.y * 80.0f - 0.0f *_integrator.y - _geomrtry_output.y- 1.0f * _phi_y + 0.0f*Aug.y);
    //_pdnn_output.z = _J_z * (-_e_R.z * 40.0f - _e_Omega.z * 80.0f - 0.0f *_integrator.z - _geomrtry_output.z - 1.0f *_phi_z + 0.0f*Aug.z); //偏航误差e_R.z很容易就趋近于0，会导致无法满足持续激励假设
    
    _pdnn_output.x = 0.01f * (-_e_R.x * 150.0f - _e_Omega.x * 40.0f - 0.0f * _integrator.x - _geomrtry_output.x - 0.0f *_phi_x + 0.0f*Aug.x); 
    _pdnn_output.y = 0.01f * (-_e_R.y * 150.0f - _e_Omega.y * 40.0f - 0.0f *_integrator.y - _geomrtry_output.y- 0.0f * _phi_y + 0.0f*Aug.y);
    _pdnn_output.z = 0.02f * (-_e_R.z * 150.0f - _e_Omega.z * 40.0f - 0.0f *_integrator.z - _geomrtry_output.z - 0.0f *_phi_z + 0.0f*Aug.z); //偏航误差e_R.z很容易就趋近于0，会导致无法满足持续激励假设
   }
    return _pdnn_output; //返回pdnn控制器输出
}

void AC_PDNN_SO3::update_i(float dt, float _ki, float _c2, float _kimax, bool limit)
{
   if (limit){

    //Vector3f delta_integrator = (_e_Omega + _e_R * _c2) * dt;
    Vector3f delta_integrator = (_e_R * _c2) * dt;
    _integrator += delta_integrator;
    
    float _integrator_x = _integrator.x;
    float _integrator_y = _integrator.y;
    float _integrator_z = _integrator.z;
    _integrator_x = constrain_float(_integrator_x, -_kimax, _kimax); //分别在xyz方向上限制积分大小
    _integrator_y = constrain_float(_integrator_y, -_kimax, _kimax);
    _integrator_z = constrain_float(_integrator_z, -_kimax, _kimax);

    _integrator.x = _integrator_x; //把限制后的值重新赋值
    _integrator.y = _integrator_y;
    _integrator.z = _integrator_z;

    _integrator =  _integrator * _ki;
} else{
    _integrator.zero();
}

}


Vector3f AC_PDNN_SO3::get_R() const
{
    return _pdnn_output_R;
}

Vector3f AC_PDNN_SO3::get_Omega() const
{
    return _pdnn_output_Omega;
}

Vector3f AC_PDNN_SO3::get_e_Omega() const
{   
    
    return _e_Omega;
}

Vector3f AC_PDNN_SO3::get_e_R() const
{   
    
    return _e_R;
}

Vector3f AC_PDNN_SO3::get_dot_Omega_c() const
{   
    
    return _dot_Omega_c;
}

Vector3f AC_PDNN_SO3::get_phi() const
{   
    Vector3f _phi;
    _phi.x = _phi_x;
    _phi.y = _phi_y;
    _phi.z = _phi_z;
    return _phi;
}

Vector3f AC_PDNN_SO3::get_J() const
{   
    Vector3f _J;
    _J.x = _J_x;
    _J.y = _J_y;
    _J.z = _J_z;
    return _J;
}

float AC_PDNN_SO3::get_Psi_R() const
{
    return _Psi_R;
}

// save_gains - save gains to eeprom
void AC_PDNN_SO3::save_gains()
{
    _kR.save();
    _kOmega.save();
    _kR_z.save();
    _kOmega_z.save();
    //_filt_E_hz.save();
    //_filt_D_hz.save();
}

// get the target filter alpha
float AC_PDNN_SO3::get_filt_E_alpha(float dt) const
{
    return calc_lowpass_alpha_dt(dt, 5.0f);
}

// get the derivative filter alpha
float AC_PDNN_SO3::get_filt_D_alpha(float dt) const
{
    return calc_lowpass_alpha_dt(dt, 5.0f);
}


  
