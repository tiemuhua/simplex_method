#include "simplex_method.h"
SimplexMethod::SolveResults SimplexMethod::solve(
    INPUT   const MatrixXd &A, const VectorXd &b_ineq, const MatrixXd &Aeq, const VectorXd &beq, const VectorXd &c, 
    OUTPUT  double &min_val, VectorXd &ori_var
) {
    init(A, b_ineq, Aeq, beq, c);
    SolveResults result = solveAssistProblem();
    if(result != SUCCESS) return result;
    solveMainProblem();
    min_val=z0_;
    return SUCCESS;
}

void SimplexMethod::init(
    INPUT   const MatrixXd &A, const VectorXd &b_ineq, const MatrixXd &Aeq, const VectorXd &beq, const VectorXd &c_ori
) {
    var_num_ = A.cols();
    ineq_num_ = b_ineq.size();
    eq_num_ = beq.size();
    row_num_=ineq_num_+eq_num_;

    A_.resize(row_num_,(var_num_<<1)+row_num_);
    b_.resize(ineq_num_ + eq_num_);
    b_.block(0, 0, ineq_num_, 1) = b_ineq;
    b_.block(ineq_num_, 0, eq_num_, 1) = beq;
    c_ = VectorXd::Zero((var_num_ << 1) + ineq_num_ + eq_num_);

    for (size_t i = 0; i < var_num_; ++i) {
        c_(i << 1) = -c_ori(i);
        c_((i << 1) + 1) = c_ori(i);
        for (size_t ineq_id = 0; ineq_id < ineq_num_; ++ineq_id) {
            A_(ineq_id,i<<1)=A(ineq_id,i);
            A_(ineq_id,(i<<1)+1)=-A(ineq_id,i);
        }
        for (size_t eq_id = 0; eq_id < eq_num_; ++eq_id) {
            A_(eq_id,i<<1)=A(eq_id,i);
            A_(eq_id,(i<<1)+1)=-A(eq_id,i);
        }
    }
    A_.block(0,var_num_<<1,row_num_,row_num_)=MatrixXd::Identity(row_num_,row_num_);

    row2var_=VectorXd::Zero(row_num_);
}

SimplexMethod::SolveResults SimplexMethod::solveAssistProblem() {
    if(eq_num_==0){
        return SUCCESS;
    }
    double min_val=0;
    VectorXd g = VectorXd::Zero((var_num_ << 1) + row_num_);
    g.block((var_num_ << 1) + ineq_num_, 0, eq_num_, 1) = -VectorXd::Ones(eq_num_);
    for (int i = ineq_num_; i < row_num_; ++i) {
        g += A_.row(i).transpose();
        min_val+=b_(i);
    }

    for(int i=0;i<row_num_;++i){
        row2var_(i)=i+(var_num_<<1);
    }
    int in = maxId(g);

    while(g(in) > 0) {
        int row = 0;
        double min_b_frac_A = inf;

        for (size_t i = 0; i < row_num_; ++i) {
            if (A_(i, in) > 0.001 && b_(i) / A_(i, in) < min_b_frac_A) {
                min_b_frac_A = b_(i) / A_(i, in);
                row = i;
            }
        }
        row2var_(row) = in;

        b_(row) /= A_(row, in);
        A_.row(row) /= A_(row, in);

        for (size_t i = 0; i < row_num_; ++i) {
            b_(i) -= A_(i, in);
            A_.row(i) -= A_(i, in) * A_.row(row);
        }

        min_val -= g(in);
        g -= g(in) * A_.row(row).transpose();

        z0_-=c_(in)*b_(row);
        c_-=c_(in)*A_.row(row).transpose();
        
        in = maxId(g);
    }

    if(min_val > 0){
        return NO_ANS;
    }

    set<int> base_ids;
    for(int i=0;i<row_num_;++i){
        base_ids.insert(row2var_(i));
    }
    set<int> erase_rows;
    int out = *(base_ids.end()--);
    while (out > (var_num_ << 1) + ineq_num_) {
        base_ids.erase(base_ids.end()--);
        int row = var2row_(out);
        auto A_row_trans = A_.row(row).transpose();
        if (A_row_trans.lpNorm<Infinity>() > 0.0001) {
            int in = maxId(A_row_trans);
            A_.row(row) /= A_(row, in);
            b_(row) /= A_(row, in);

            for (size_t i = 0; i < row_num_; ++i) {
                A_.row(i) -= A_(i, in) * A_.row(row);
                b_(i) -= A_(i, in);
            }

            c_ -= c_(in) * A_.row(row).transpose();
            z0_ -= c_(in) * b_(row);
        } else {
            erase_rows.insert(row);
        }
        out = *(base_ids.end()--);
    }
    int new_row = 0;
    MatrixXd new_mat = A_;
    VectorXd new_row2var = row2var_;
    VectorXd new_b = b_;
    for (int i = 0; i < A_.rows(); ++i) {
        if (!erase_rows.count(i)) {
            new_mat.row(new_row) = A_.row(i);
            new_row2var(new_row) = row2var_(i);
            new_b(new_row) = b_(i);
            new_row++;
        }else{
            erase_rows.erase(i);
        }
    }
    A_ = new_mat.block(0, 0, new_row, (var_num_ << 1) + ineq_num_);
    b_ = new_b.block(0, 0, new_row, 1);
    row2var_ = new_row2var.block(0, 0, new_row, 1);

    return SUCCESS;
}

void SimplexMethod::solveMainProblem() {
    int in = maxId(c_);
    cout<<"in:\t"<<in<<"c_.size"<<c_.size()<<endl;
    while(c_(in) > 0){
        int row = 0;
        double min_b_frac_A = inf;
        for (size_t i = 0; i < row_num_; ++i) {
            if (A_(i, in) > 0.001 && b_(i) / A_(i, in) < min_b_frac_A) {
                min_b_frac_A = b_(i) / A_(i, in);
                row = i;
            }
        }
        row2var_(row) = in;

        b_(row) /= A_(row, in);
        A_.row(row) /= A_(row, in);

        for (size_t i = 0; i < row_num_; ++i) {
            if(i==row)continue;
            b_(i) -= A_(i, in)*b_(row);
            A_.row(i) -= A_(i, in) * A_.row(row);
        }

        z0_ -= c_(in) * b_(row);
        c_ -= c_(in) * A_.row(row).transpose();
        cout<<c_.transpose()<<"\n\n";
        cout << z0_ << endl;
        in = maxId(c_);
    }
}
inline int SimplexMethod::maxId(const VectorXd &vec){
    double max_val=-inf;
    int max_id=0;
    for(int i=0;i<vec.size();++i){
        if(max_val<vec(i)){
            max_val=vec(i);
            max_id=i;
        }
    }
    return max_id;
}
