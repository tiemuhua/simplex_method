#include <bits/stdc++.h>

#include <Eigen/Core>
using namespace Eigen;
using namespace std;

#define INPUT
#define OUTPUT

class SimplexMethod {
   public:
    enum SolveResults { SUCCESS, NO_ANS, INF };

   private:
    const double inf = 0x3f3f3f3f;
    set<int> baseVars_;
    MatrixXd A_;
    VectorXd b_;
    VectorXd c_;
    double z0_;

    VectorXd row2var_;
    VectorXd var2row_;
    MatrixXd x2ori_var;

    int var_num_;
    int ineq_num_;
    int eq_num_;
    int row_num_;

    void clear();
    void init(
        INPUT   const MatrixXd &A_ineq, const VectorXd &b_ineq, const MatrixXd &Aeq, const VectorXd &beq, const VectorXd &c_ori
    );
    SolveResults solveAssistProblem();
    void solveMainProblem();
    inline int maxId(const VectorXd &vec);
    inline int min_BifracAik_id(const int in_id);

   public:
    SolveResults solve(
        INPUT   const MatrixXd &A, const VectorXd &b, const MatrixXd &Aeq, const VectorXd &beq, const VectorXd &c,
        OUTPUT  double &min_val, VectorXd &ori_var
    );
    SolveResults solve(
        INPUT   const MatrixXd &A, const VectorXd &b, const VectorXd &c,
        OUTPUT  double &min_val, VectorXd &ori_var
    );
    SimplexMethod(){}
    ~SimplexMethod(){}
};