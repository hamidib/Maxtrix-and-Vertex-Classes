//============================================================
// file: main.cpp
// HW by Bijan A. Hamidi
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>
#include "matrix_3dT.h"

#include "vector_3dT.h"
template <typename T>
void print(T v) {
    std::cout << v << std::endl;
}
template <typename T>
void show_vect(T v) {
    std::cout << v.name() << " is: " << v << std::endl;
}
template <typename T>
void show_mat(T m) {
    std::cout << m.name() << " is: " << m << std::endl;
}
void test_vectors() {
    print("\n====================  TESTING VECTORS  ========================");
    vector3dD u("u", 3, {1,  2,  4});
    vector3dD v("v", 3, {8, 16, 32});
    vector3dD i("i", 3, {1, 0, 0}), j("j", 3, {0, 1, 0}), k("k", 3, {0, 0, 1});
    vector3dD w(3 * i + 4 * j - 2 * k);
    show_vect(u);
    show_vect(v);
    show_vect(i);
    show_vect(j);
    show_vect(k);
    show_vect(w);
    assert(u == u);
    assert(u != v);
    assert(u + v == v + u);
    assert(u - v == -(v - u));
    assert(-(-u) == u);
    assert(3.0 + u == u + 3.0);
    assert(3.0 * u == u * 3.0);
    assert((u - 3.0) == -(3.0 - u));
    assert((5.0 * u) / 5.0 == u);
    assert(u + vector3dD::zero() == u);
    assert(i.dot(j) == j.dot(k) == k.dot(i) == 0);
    assert(i.cross(j) == k);
    assert(j.cross(k) == i);
    assert(k.cross(i) == j);
    assert(u.cross(v) == -v.cross(u));
    assert(u.cross(v + w) == u.cross(v) + u.cross(w));
    assert((u.cross(v)).dot(u) == 0);
    print(i.angle(j));
    print(M_PI/2);
    assert(i.angle(j) == M_PI_2); //means pi over 2
    assert(j.angle(k) == M_PI_2);
    assert(k.angle(i) == M_PI_2);
    vector3dD uhat = u / u.magnitude();
    show_vect(u);
    show_vect(uhat);
    print(uhat.magnitude());
    assert(uhat.magnitude() - 1.0 < 1.0e-10);
    print("...test vectors assertions passed");
    print("====================  FINISHED testing vectors  ========================");
}

void test_matrices() {
    print("\n====================  TESTING MATRICES  ========================");
    matrix3dD a("a", 3, {3, 2, 0, 0, 0, 1, 2, -2, 1});
    matrix3dD b("b", 3, {1, 0, 5, 2, 1, 6, 3,  4, 0});
    matrix3dD ainv = a.inverse();//.transpose();
    matrix3dD binv = b.inverse();
    print(a);
    print(ainv);
    print(a * ainv);
    
//    print(b);
//    print(a*b);
//
//    print(ainv);
    print(binv);
//    //print( b + a);//test
//    //print(a - a); //test
//    //print(a / b); //test
//    //print(a * b);//test
    matrix3dD btimes3 = b * 3;
//    print(b);
//    print(btimes3);
    matrix3dD btimesa = b * a;
//    print(btimesa);
    matrix3dD atimesb = a * b;
//    print(atimesb);
//    print(a * ainv);
//    print(b * binv);
    assert(a * ainv == matrix3dD::identity(3)); //not exactly equal to zero//
//    print(a.cofactor());
//    print(b.cofactor());
//    print(a.adjugate());
//    print(b.adjugate());
//    print(a.transpose());
//    print(b.transpose());
//    print(a);
////    print(b);
//    std::cout << "==================\n";
//    print(ainv.transpose());
////    print(binv);
//    std::cout << "==================\n";
//    print(a*(ainv.transpose()));
////    print(ainv*a);
    assert(a * ainv == ainv * a);//
    assert(b * binv == matrix3dD::identity(3));
    assert(b * binv == binv * b);
    assert(a.transpose().transpose() == a); //transpose of transpose is normal
    assert(a.transpose().determinant() == a.determinant());
    assert(a + b == b + a);
    assert(a - b == -(b - a));
    assert(3.0 + a == a + 3.0);
    assert(3.0 * a == a * 3.0);
    assert((a + 3.0) - 3.0 == a);
    assert((3.0 * a) / 3.0 == a);
    assert(-(-a) == a);
    matrix3dD zerod("zerod", 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    assert(zerod.determinant() == 0);
    print("...test matrices assertions passed");
    print("====================  FINISHED testing matrices  ========================");
}
void test_matrices_and_vectors() {
    print("\n====================  TESTING MATRICES and VECTORS  ========================");
    vector3dD p("p", 2, {1, 2});
    matrix3dD m("m", 2, {1, 2, 3, 4});
    show_vect(p);
    show_mat(m);
    //assert(p * m == m * p);
    vector3dD q("q", 3, {1, 2, 3});
    matrix3dD n("n", 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    show_vect(q);
    show_mat(n);
    //assert(q * n == n * q);
    print("...test_matrices_and_vectors assertions passed");
    print("====================  FINISHED testing matrices and vectors  ========================");
}
int main(int argc, const char * argv[]) {
    test_vectors();
    test_matrices();
    test_matrices_and_vectors();
    print("... program completed...\n");
    return 0;
}

/* SAMPLE OUTPUT
 ====================  TESTING VECTORS  ========================
 u is: <'u', 1 2 4 0>
 v is: <'v', 8 16 32 0>
 i is: <'i', 1 0 0 0>
 j is: <'j', 0 1 0 0>
 k is: <'k', 0 0 1 0>
 3.000000*i+4.000000*j-2.000000*k is: <'3.000000*i+4.000000*j-2.000000*k', 3 4 -2 0> 1.5708
 1.5708
 u is: <'u', 1 2 4 0>
 u/4.582576 is: <'u/4.582576', 0.218218 0.436436 0.872872 0>
 1
 ...test vectors assertions passed
 ==================== FINISHED testing vectors ========================
 ====================  TESTING MATRICES  ========================
 <'a', <'col0', 3 2 0 0><'col1', 0 0 1 0><'col2', 2 -2 1 0>> OR by rows... 
 3 02
 2 0 -2
 0 11
 >
 <'b', <'col0', 1 0 5 0><'col1', 2 1 6 0><'col2', 3 4 0 0>> OR by rows...
 1 23
 0 14
 5 60
 >
 <'Co(a)T/10.000000', <'col0/10.000000', 0.2 -0.2 0.2 0><'col1/10.000000', 0.2 0.3 -0.3 0><'col2/10.000000', -0 1 0 0>> OR by rows...
 0.2 0.2 -0
 -0.2 0.3 1
 0.2 -0.3 0
 >
 <'Co(b)T/1.000000', <'col0/1.000000', -24 20 -5 0><'col1/1.000000', 18 -15 4 0><'col2/1.000000', 5 -4 1 0>> OR by rows...
 -24 18 5
 20 -15 -4
 -5 4 1
 >
 <'a*Co(a)T/10.000000', <'col0', 1 0 0 0><'col1', 1.11022e-16 1 0 0><'col2', 0 0 1 0>> OR by rows...
 1 1.11022e-16 0
 0 10
 0 01
 >
 <'b*Co(b)T/1.000000', <'col0', 1 0 0 0><'col1', 0 1 0 0><'col2', 0 0 1 0>> OR by rows...
 1 00
 0 10
 0 01
 >
 ...test matrices assertions passed
 ==================== FINISHED testing matrices ========================
 ====================  TESTING MATRICES and VECTORS  ========================
 p is: <'p', 1 2 0>
 m is: <'m', <'col0', 1 2 3 0><'col1', 4 0 0 0><'col2', 0 0 0 0>> OR by rows...
 1 40
 2 00
 3 00
 >
 q is: <'q', 1 2 3 0>
 n is: <'n', <'col0', 1 2 3 0><'col1', 4 5 6 0><'col2', 7 8 9 0>> OR by rows...
 1 47
 2 58
 3 69
 >
 ...test_matrices_and_vectors assertions passed
 ==================== FINISHED testing matrices and vectors ======================== ... program completed...
 Program ended with exit code: 0
*/