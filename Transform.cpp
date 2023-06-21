#include "Transform.h"
  
mat3 Transform::rotate(const float degrees, const vec3& axis)
{
    const double rads = glm::radians(degrees);
    const float sinTheta = sin(rads);
    const float cosTheta = cos(rads);
    const float x = axis.x;
    const float y = axis.y;
    const float z = axis.z;
    return cosTheta * glm::mat3()
        + (1 - cosTheta) * glm::mat3(x * x, x * y, x * z, y * x, y * y, y * z, z * x, z * y, z * z)
        + sinTheta * glm::mat3(0, z, -y, -z, 0, x, y, -x, 0);
}

void Transform::left(float degrees, vec3& eye, vec3& up)
{
    mat3 rot = rotate(degrees, up);
    eye = rot * eye;
}

void Transform::up(float degrees, vec3& eye, vec3& up)
{
    vec3 axis = normalize(cross(eye, up));
    mat3 rot = rotate(degrees, axis);
    eye = rot * eye;
    up = rot * up;
}

mat4 Transform::lookAt(const vec3& eye, const vec3& center, const vec3& up)
{
    vec3 w = normalize(eye);
    vec3 u = normalize(cross(up, w));
    vec3 v = cross(w, u);

    mat4 m1(u.x, v.x, w.x, 0, u.y, v.y, w.y, 0, u.z, v.z, w.z, 0, 0, 0, 0, 1);
    mat4 m2(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -eye.x, -eye.y, -eye.z, 1);
    return m1 * m2;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{
    float theta = fovy * pi / 180 / 2;
    float d = cos(theta) / sin(theta);
    float A = -(zFar + zNear) / (zFar - zNear);
    float B = -2 * (zFar * zNear) / (zFar - zNear);
    mat4 ret(d / aspect, 0, 0, 0,
        0, d, 0, 0,
        0, 0, A, -1,
        0, 0, B, 0);
    return ret;
}

mat4 Transform::scale(const float& sx, const float& sy, const float& sz)
{
    mat4 ret(sx, 0, 0, 0,
        0, sy, 0, 0,
        0, 0, sz, 0,
        0, 0, 0, 1);
    return ret;
}

mat4 Transform::translate(const float& tx, const float& ty, const float& tz)
{
    mat4 ret(1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        tx, ty, tz, 1);
    return ret;
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3& up, const vec3& zvec)
{
    vec3 x = glm::cross(up, zvec);
    vec3 y = glm::cross(zvec, x);
    vec3 ret = glm::normalize(y);
    return ret;
}


Transform::Transform()
{

}

Transform::~Transform()
{

}
