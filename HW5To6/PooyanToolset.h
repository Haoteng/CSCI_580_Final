//
//  PooyanToolset.h
//  CSCI_580_HW1
//
//  Created by Pooyan Behnamghader on 2/9/15.
//  Copyright (c) 2015 Pooyan Behnamghader. All rights reserved.
//

#ifndef CSCI_580_HW1_PooyanToolset_h
#define CSCI_580_HW1_PooyanToolset_h

#define MIN(A,B) ((A)>(B)?(B):(A))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A) ((A)>0?(1):(-1))
#define ABS(A) ((A)>0?(A):(-1*(A)))

#ifndef MAX_INT
#define MAX_INT 2147483647
#endif


class PooyanToolset {
    public:
    // This function checks if a number is between two numbers
    template <typename Type>
    static bool isInRange(Type, Type, Type);
    // This function fit a value in a min max frame
    template <typename Type>
    static void fitValue (Type&, Type, Type);
    // This is copying from a void array to a type
    template <typename Type>
    static void copyFromVoid(Type&, void*);
    // This function finds A, B, C in Ax+By+C=0
    static void lineEquation (float [3],float [2] , float[2]);
    // This function finds A, B, C, D in Ax+By+Cz+D=0 ;
    static void plateEquation (float [4],float[3],float[3],float[3]);
    // This function finds z giving x,y in Ax+By+Cz+D=0
    static float calZ (float [4],float,float);
    // This function convert degree to radian
    static float degreeToRadian(float);
    // This function calculates the length of a vector
    static float lengthCalculation (float [3]) ;
    // This function normalize a vector
    static int normalize (float [3]) ;
    // This function check the equlity of two vectors
    static bool isEqual (float [3], float[3]);
    // This function calculates inner production of two vectors
    static float innerProduct (float [3],float[3]) ;
    // This function cross products two matrices;
    static void crossProduct (float [3],float[3],float[3]) ;
    // This function is  end*t+(1-t)*start;
    static float relativeNum(float,float,float);
    // This function calculate bilinear interpolation considering u and v in [0,1];
    static float bilinear (float,float,float,float,float,float);

};

template <typename Type>
bool PooyanToolset::isInRange(Type x, Type min, Type max){
    return x>=min && x<=max ;
}

template <typename Type>
void PooyanToolset::fitValue (Type &x, Type min, Type max){
    if (x>max)
        x=max;
    if (x<min)
        x=min;
}

template <typename Type>
void PooyanToolset::copyFromVoid(Type& d, void *p) {
	for (int i=0;i<sizeof(Type);i++)
		((char *)&d)[i]=((char *)p)[i];
}

#endif
