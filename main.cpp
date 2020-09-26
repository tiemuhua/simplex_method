#include <bits/stdc++.h>

#include <Eigen/Core>

#include "simplex_method.h"

using namespace std;
using namespace Eigen;
int main() {
    freopen("../in", "r", stdin);
    VectorXd c, b, beq;
    MatrixXd A, Aeq;
    int var_num, ineq_num, eq_num;
    double tmp;
    cin >> var_num >> ineq_num >> eq_num;
    A.resize(ineq_num,var_num);
    b.resize(ineq_num);
    Aeq.resize(eq_num,var_num);
    beq.resize(eq_num);
    c.resize(var_num);

    for(int i=0;i<ineq_num;++i){
        for(int j=0;j<var_num;++j){
            cin>>tmp;
            A(i,j)=tmp;
        }
    }
    for(int i=0;i<ineq_num;++i){
        cin>>tmp;
        b(i)=tmp;
    }
    for(int i=0;i<eq_num;++i){
        for(int j=0;j<var_num;++j){
            cin>>tmp;
            Aeq(i,j)=tmp;
        }
    }
    for(int i=0;i<eq_num;++i){
        cin>>tmp;
        beq(i)=tmp;
    }
    for(int i=0;i<var_num;++i){
        cin>>tmp;
        c(i)=tmp;
    }

    SimplexMethod simplex_method;
    double min_val = 0;
    VectorXd ans_vec;
    SimplexMethod::SolveResults result = simplex_method.solve(A, b, Aeq, beq, c, min_val, ans_vec);
    if (result == SimplexMethod::NO_ANS) {
        cout << "no ans" << endl;
    }
    cout << "min val:\t" << min_val << endl;
    cout<<ans_vec.transpose()<<endl;
}