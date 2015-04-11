#include "PooyanToolset.h"
#include <math.h>

void PooyanToolset::lineEquation(float coef[3], float p1[2], float p2[2]) {
	coef[0] = p1[1] - p2[1];
	coef[1] = p2[0] - p1[0];
	coef[2] = -coef[1] * p1[1] - coef[0] * p1[0];
}

void PooyanToolset::plateEquation(float coef[4], float p1[3], float p2[3],
		float p3[3]) {
	coef[0] = p2[1] * (p1[2] - p3[2]) + p1[1] * (p3[2] - p2[2])
			+ p3[1] * (p2[2] - p1[2]);
	coef[1] = p2[2] * (p1[0] - p3[0]) + p1[2] * (p3[0] - p2[0])
			+ p3[2] * (p2[0] - p1[0]);
	coef[2] = p2[0] * (p1[1] - p3[1]) + p1[0] * (p3[1] - p2[1])
			+ p3[0] * (p2[1] - p1[1]);
	coef[3] = -coef[0] * p1[0] - coef[1] * p1[1] - coef[2] * p1[2];
	/*
	float min = MIN(ABS(coef[0]),MIN(ABS(coef[1]),MIN(ABS(coef[2]),ABS(coef[3]))));
	min=MAX(min,0.00001);
	for (int i=0 ; i<4;i++)
		coef[i]/=min;
	*/
}

float PooyanToolset::degreeToRadian(float degree) {
	return degree * M_PI / 180.0;
}

float PooyanToolset::lengthCalculation(float v[3]) {
	return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
}

int PooyanToolset::normalize(float v[3]) {
	float length = lengthCalculation(v);
	if (length == 0)
		return 1;
	v[0] /= length;
	v[1] /= length;
	v[2] /= length;
	return 0;
}

bool PooyanToolset::isEqual(float v1[3], float v2[3]) {
	return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

float PooyanToolset::innerProduct(float v1[3], float v2[3]) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void PooyanToolset::crossProduct(float v1[3], float v2[3], float r[3]) {
	r[0] = v1[1] * v2[2] - v1[2] * v2[1];
	r[1] = v1[2] * v2[0] - v1[0] * v2[2];
	r[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
float PooyanToolset::relativeNum(float start, float end, float t) {
	return end * t + (1 - t) * start;
}

float PooyanToolset::calZ(float coef[4], float x, float y) {
	return (-coef[0] * x - coef[1] * y - coef[3]) / coef[2];
}

float PooyanToolset::bilinear(float u,float v,float V11,float V21, float V22, float V12){
	return V11*(1-u)*(1-v)+V21*u*(1-v)+V22*u*v+V12*(1-u)*v;
}
