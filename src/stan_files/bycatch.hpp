/*
    bycatch is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bycatch is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bycatch.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.0

#include <stan/model/model_header.hpp>

namespace model_bycatch_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_bycatch");
    reader.add_event(74, 72, "end", "model_bycatch");
    return reader;
}

#include <meta_header.hpp>
 class model_bycatch : public prob_grad {
private:
    int n_row;
    vector_d effort;
    vector<int> yint;
    vector<int> time;
    int n_year;
    int K;
    matrix_d x;
    int family;
    int time_varying;
public:
    model_bycatch(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_bycatch(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_bycatch_namespace::model_bycatch";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "n_row", "int", context__.to_vec());
            n_row = int(0);
            vals_i__ = context__.vals_i("n_row");
            pos__ = 0;
            n_row = vals_i__[pos__++];
            current_statement_begin__ = 3;
            validate_non_negative_index("effort", "n_row", n_row);
            context__.validate_dims("data initialization", "effort", "vector_d", context__.to_vec(n_row));
            validate_non_negative_index("effort", "n_row", n_row);
            effort = vector_d(static_cast<Eigen::VectorXd::Index>(n_row));
            vals_r__ = context__.vals_r("effort");
            pos__ = 0;
            size_t effort_i_vec_lim__ = n_row;
            for (size_t i_vec__ = 0; i_vec__ < effort_i_vec_lim__; ++i_vec__) {
                effort[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 4;
            validate_non_negative_index("yint", "n_row", n_row);
            context__.validate_dims("data initialization", "yint", "int", context__.to_vec(n_row));
            validate_non_negative_index("yint", "n_row", n_row);
            yint = std::vector<int>(n_row,int(0));
            vals_i__ = context__.vals_i("yint");
            pos__ = 0;
            size_t yint_limit_0__ = n_row;
            for (size_t i_0__ = 0; i_0__ < yint_limit_0__; ++i_0__) {
                yint[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("time", "n_row", n_row);
            context__.validate_dims("data initialization", "time", "int", context__.to_vec(n_row));
            validate_non_negative_index("time", "n_row", n_row);
            time = std::vector<int>(n_row,int(0));
            vals_i__ = context__.vals_i("time");
            pos__ = 0;
            size_t time_limit_0__ = n_row;
            for (size_t i_0__ = 0; i_0__ < time_limit_0__; ++i_0__) {
                time[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "n_year", "int", context__.to_vec());
            n_year = int(0);
            vals_i__ = context__.vals_i("n_year");
            pos__ = 0;
            n_year = vals_i__[pos__++];
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            current_statement_begin__ = 8;
            validate_non_negative_index("x", "n_row", n_row);
            validate_non_negative_index("x", "K", K);
            context__.validate_dims("data initialization", "x", "matrix_d", context__.to_vec(n_row,K));
            validate_non_negative_index("x", "n_row", n_row);
            validate_non_negative_index("x", "K", K);
            x = matrix_d(static_cast<Eigen::VectorXd::Index>(n_row),static_cast<Eigen::VectorXd::Index>(K));
            vals_r__ = context__.vals_r("x");
            pos__ = 0;
            size_t x_m_mat_lim__ = n_row;
            size_t x_n_mat_lim__ = K;
            for (size_t n_mat__ = 0; n_mat__ < x_n_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < x_m_mat_lim__; ++m_mat__) {
                    x(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "family", "int", context__.to_vec());
            family = int(0);
            vals_i__ = context__.vals_i("family");
            pos__ = 0;
            family = vals_i__[pos__++];
            current_statement_begin__ = 10;
            context__.validate_dims("data initialization", "time_varying", "int", context__.to_vec());
            time_varying = int(0);
            vals_i__ = context__.vals_i("time_varying");
            pos__ = 0;
            time_varying = vals_i__[pos__++];

            // validate, data variables
            current_statement_begin__ = 2;
            check_greater_or_equal(function__,"n_row",n_row,0);
            current_statement_begin__ = 3;
            current_statement_begin__ = 4;
            current_statement_begin__ = 5;
            current_statement_begin__ = 6;
            check_greater_or_equal(function__,"n_year",n_year,0);
            current_statement_begin__ = 7;
            check_greater_or_equal(function__,"K",K,0);
            current_statement_begin__ = 8;
            current_statement_begin__ = 9;
            current_statement_begin__ = 10;
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 13;
            validate_non_negative_index("beta", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 14;
            validate_non_negative_index("est_time_dev", "(time_varying * (n_year - 1))", (time_varying * (n_year - 1)));
            num_params_r__ += (time_varying * (n_year - 1));
            current_statement_begin__ = 15;
            validate_non_negative_index("sigma_rw", "time_varying", time_varying);
            num_params_r__ += time_varying;
            current_statement_begin__ = 16;
            validate_non_negative_index("nb2_phi", "family", family);
            num_params_r__ += family;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_bycatch() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("beta")))
            throw std::runtime_error("variable beta missing");
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "K", K);
        context__.validate_dims("initialization", "beta", "vector_d", context__.to_vec(K));
        vector_d beta(static_cast<Eigen::VectorXd::Index>(K));
        for (int j1__ = 0U; j1__ < K; ++j1__)
            beta(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable beta: ") + e.what());
        }

        if (!(context__.contains_r("est_time_dev")))
            throw std::runtime_error("variable est_time_dev missing");
        vals_r__ = context__.vals_r("est_time_dev");
        pos__ = 0U;
        validate_non_negative_index("est_time_dev", "(time_varying * (n_year - 1))", (time_varying * (n_year - 1)));
        context__.validate_dims("initialization", "est_time_dev", "vector_d", context__.to_vec((time_varying * (n_year - 1))));
        vector_d est_time_dev(static_cast<Eigen::VectorXd::Index>((time_varying * (n_year - 1))));
        for (int j1__ = 0U; j1__ < (time_varying * (n_year - 1)); ++j1__)
            est_time_dev(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(est_time_dev);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable est_time_dev: ") + e.what());
        }

        if (!(context__.contains_r("sigma_rw")))
            throw std::runtime_error("variable sigma_rw missing");
        vals_r__ = context__.vals_r("sigma_rw");
        pos__ = 0U;
        validate_non_negative_index("sigma_rw", "time_varying", time_varying);
        context__.validate_dims("initialization", "sigma_rw", "double", context__.to_vec(time_varying));
        std::vector<double> sigma_rw(time_varying,double(0));
        for (int i0__ = 0U; i0__ < time_varying; ++i0__)
            sigma_rw[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < time_varying; ++i0__)
            try {
            writer__.scalar_lb_unconstrain(0,sigma_rw[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_rw: ") + e.what());
        }

        if (!(context__.contains_r("nb2_phi")))
            throw std::runtime_error("variable nb2_phi missing");
        vals_r__ = context__.vals_r("nb2_phi");
        pos__ = 0U;
        validate_non_negative_index("nb2_phi", "family", family);
        context__.validate_dims("initialization", "nb2_phi", "double", context__.to_vec(family));
        std::vector<double> nb2_phi(family,double(0));
        for (int i0__ = 0U; i0__ < family; ++i0__)
            nb2_phi[i0__] = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < family; ++i0__)
            try {
            writer__.scalar_lb_unconstrain(0,nb2_phi[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable nb2_phi: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(K,lp__);
            else
                beta = in__.vector_constrain(K);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  est_time_dev;
            (void) est_time_dev;  // dummy to suppress unused var warning
            if (jacobian__)
                est_time_dev = in__.vector_constrain((time_varying * (n_year - 1)),lp__);
            else
                est_time_dev = in__.vector_constrain((time_varying * (n_year - 1)));

            vector<local_scalar_t__> sigma_rw;
            size_t dim_sigma_rw_0__ = time_varying;
            sigma_rw.reserve(dim_sigma_rw_0__);
            for (size_t k_0__ = 0; k_0__ < dim_sigma_rw_0__; ++k_0__) {
                if (jacobian__)
                    sigma_rw.push_back(in__.scalar_lb_constrain(0,lp__));
                else
                    sigma_rw.push_back(in__.scalar_lb_constrain(0));
            }

            vector<local_scalar_t__> nb2_phi;
            size_t dim_nb2_phi_0__ = family;
            nb2_phi.reserve(dim_nb2_phi_0__);
            for (size_t k_0__ = 0; k_0__ < dim_nb2_phi_0__; ++k_0__) {
                if (jacobian__)
                    nb2_phi.push_back(in__.scalar_lb_constrain(0,lp__));
                else
                    nb2_phi.push_back(in__.scalar_lb_constrain(0));
            }


            // transformed parameters
            current_statement_begin__ = 19;
            validate_non_negative_index("lambda", "n_row", n_row);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda(static_cast<Eigen::VectorXd::Index>(n_row));
            (void) lambda;  // dummy to suppress unused var warning

            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda,DUMMY_VAR__);
            current_statement_begin__ = 20;
            validate_non_negative_index("pred", "n_row", n_row);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  pred(static_cast<Eigen::VectorXd::Index>(n_row));
            (void) pred;  // dummy to suppress unused var warning

            stan::math::initialize(pred, DUMMY_VAR__);
            stan::math::fill(pred,DUMMY_VAR__);
            current_statement_begin__ = 21;
            validate_non_negative_index("time_dev", "(time_varying * n_year)", (time_varying * n_year));
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  time_dev(static_cast<Eigen::VectorXd::Index>((time_varying * n_year)));
            (void) time_dev;  // dummy to suppress unused var warning

            stan::math::initialize(time_dev, DUMMY_VAR__);
            stan::math::fill(time_dev,DUMMY_VAR__);


            current_statement_begin__ = 22;
            stan::math::assign(pred, multiply(x,beta));
            current_statement_begin__ = 24;
            if (as_bool(logical_eq(time_varying,1))) {

                current_statement_begin__ = 25;
                stan::model::assign(time_dev, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            0, 
                            "assigning variable time_dev");
                current_statement_begin__ = 26;
                for (int i = 2; i <= n_year; ++i) {

                    current_statement_begin__ = 27;
                    stan::model::assign(time_dev, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                get_base1(est_time_dev,(i - 1),"est_time_dev",1), 
                                "assigning variable time_dev");
                }
            }
            current_statement_begin__ = 31;
            for (int i = 1; i <= n_row; ++i) {

                current_statement_begin__ = 32;
                if (as_bool(logical_eq(time_varying,1))) {

                    current_statement_begin__ = 33;
                    stan::model::assign(pred, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::model::deep_copy((get_base1(pred,i,"pred",1) + (time_varying * get_base1(time_dev,get_base1(time,i,"time",1),"time_dev",1)))), 
                                "assigning variable pred");
                }
                current_statement_begin__ = 35;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::exp((get_base1(pred,i,"pred",1) + stan::math::log(get_base1(effort,i,"effort",1)))), 
                            "assigning variable lambda");
            }

            // validate transformed parameters
            for (int i0__ = 0; i0__ < n_row; ++i0__) {
                if (stan::math::is_uninitialized(lambda(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: lambda" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < n_row; ++i0__) {
                if (stan::math::is_uninitialized(pred(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: pred" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < (time_varying * n_year); ++i0__) {
                if (stan::math::is_uninitialized(time_dev(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: time_dev" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 19;
            current_statement_begin__ = 20;
            current_statement_begin__ = 21;

            // model body

            current_statement_begin__ = 39;
            lp_accum__.add(student_t_log<propto__>(beta, 3, 0, 2));
            current_statement_begin__ = 41;
            if (as_bool(logical_eq(time_varying,1))) {

                current_statement_begin__ = 42;
                lp_accum__.add(student_t_log<propto__>(sigma_rw, 3, 0, 1));
                current_statement_begin__ = 43;
                lp_accum__.add(student_t_log<propto__>(get_base1(est_time_dev,1,"est_time_dev",1), 3, 0, 2));
                current_statement_begin__ = 44;
                for (int i = 2; i <= (n_year - 1); ++i) {

                    current_statement_begin__ = 46;
                    lp_accum__.add(normal_log<propto__>(get_base1(est_time_dev,i,"est_time_dev",1), get_base1(est_time_dev,(i - 1),"est_time_dev",1), get_base1(sigma_rw,1,"sigma_rw",1)));
                }
            }
            current_statement_begin__ = 50;
            if (as_bool(logical_eq(family,0))) {

                current_statement_begin__ = 51;
                for (int i = 1; i <= n_row; ++i) {

                    current_statement_begin__ = 52;
                    lp_accum__.add(poisson_log<propto__>(get_base1(yint,i,"yint",1), get_base1(lambda,i,"lambda",1)));
                }
            } else {

                current_statement_begin__ = 55;
                lp_accum__.add(student_t_log<propto__>(nb2_phi, 3, 0, 2));
                current_statement_begin__ = 56;
                for (int i = 1; i <= n_row; ++i) {

                    current_statement_begin__ = 57;
                    lp_accum__.add(neg_binomial_2_log<propto__>(get_base1(yint,i,"yint",1), get_base1(lambda,i,"lambda",1), get_base1(nb2_phi,1,"nb2_phi",1)));
                }
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta");
        names__.push_back("est_time_dev");
        names__.push_back("sigma_rw");
        names__.push_back("nb2_phi");
        names__.push_back("lambda");
        names__.push_back("pred");
        names__.push_back("time_dev");
        names__.push_back("log_lik");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((time_varying * (n_year - 1)));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(time_varying);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(family);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_row);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_row);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((time_varying * n_year));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_row);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_bycatch_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d beta = in__.vector_constrain(K);
        vector_d est_time_dev = in__.vector_constrain((time_varying * (n_year - 1)));
        vector<double> sigma_rw;
        size_t dim_sigma_rw_0__ = time_varying;
        for (size_t k_0__ = 0; k_0__ < dim_sigma_rw_0__; ++k_0__) {
            sigma_rw.push_back(in__.scalar_lb_constrain(0));
        }
        vector<double> nb2_phi;
        size_t dim_nb2_phi_0__ = family;
        for (size_t k_0__ = 0; k_0__ < dim_nb2_phi_0__; ++k_0__) {
            nb2_phi.push_back(in__.scalar_lb_constrain(0));
        }
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(beta[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < (time_varying * (n_year - 1)); ++k_0__) {
            vars__.push_back(est_time_dev[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < time_varying; ++k_0__) {
            vars__.push_back(sigma_rw[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < family; ++k_0__) {
            vars__.push_back(nb2_phi[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            current_statement_begin__ = 19;
            validate_non_negative_index("lambda", "n_row", n_row);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda(static_cast<Eigen::VectorXd::Index>(n_row));
            (void) lambda;  // dummy to suppress unused var warning

            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda,DUMMY_VAR__);
            current_statement_begin__ = 20;
            validate_non_negative_index("pred", "n_row", n_row);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  pred(static_cast<Eigen::VectorXd::Index>(n_row));
            (void) pred;  // dummy to suppress unused var warning

            stan::math::initialize(pred, DUMMY_VAR__);
            stan::math::fill(pred,DUMMY_VAR__);
            current_statement_begin__ = 21;
            validate_non_negative_index("time_dev", "(time_varying * n_year)", (time_varying * n_year));
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  time_dev(static_cast<Eigen::VectorXd::Index>((time_varying * n_year)));
            (void) time_dev;  // dummy to suppress unused var warning

            stan::math::initialize(time_dev, DUMMY_VAR__);
            stan::math::fill(time_dev,DUMMY_VAR__);


            current_statement_begin__ = 22;
            stan::math::assign(pred, multiply(x,beta));
            current_statement_begin__ = 24;
            if (as_bool(logical_eq(time_varying,1))) {

                current_statement_begin__ = 25;
                stan::model::assign(time_dev, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            0, 
                            "assigning variable time_dev");
                current_statement_begin__ = 26;
                for (int i = 2; i <= n_year; ++i) {

                    current_statement_begin__ = 27;
                    stan::model::assign(time_dev, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                get_base1(est_time_dev,(i - 1),"est_time_dev",1), 
                                "assigning variable time_dev");
                }
            }
            current_statement_begin__ = 31;
            for (int i = 1; i <= n_row; ++i) {

                current_statement_begin__ = 32;
                if (as_bool(logical_eq(time_varying,1))) {

                    current_statement_begin__ = 33;
                    stan::model::assign(pred, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::model::deep_copy((get_base1(pred,i,"pred",1) + (time_varying * get_base1(time_dev,get_base1(time,i,"time",1),"time_dev",1)))), 
                                "assigning variable pred");
                }
                current_statement_begin__ = 35;
                stan::model::assign(lambda, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::exp((get_base1(pred,i,"pred",1) + stan::math::log(get_base1(effort,i,"effort",1)))), 
                            "assigning variable lambda");
            }

            // validate transformed parameters
            current_statement_begin__ = 19;
            current_statement_begin__ = 20;
            current_statement_begin__ = 21;

            // write transformed parameters
            if (include_tparams__) {
            for (int k_0__ = 0; k_0__ < n_row; ++k_0__) {
            vars__.push_back(lambda[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < n_row; ++k_0__) {
            vars__.push_back(pred[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < (time_varying * n_year); ++k_0__) {
            vars__.push_back(time_dev[k_0__]);
            }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 62;
            validate_non_negative_index("log_lik", "n_row", n_row);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  log_lik(static_cast<Eigen::VectorXd::Index>(n_row));
            (void) log_lik;  // dummy to suppress unused var warning

            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik,DUMMY_VAR__);


            current_statement_begin__ = 63;
            if (as_bool(logical_eq(family,0))) {

                current_statement_begin__ = 64;
                for (int n = 1; n <= n_row; ++n) {

                    current_statement_begin__ = 65;
                    stan::model::assign(log_lik, 
                                stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                                poisson_log(get_base1(yint,n,"yint",1),get_base1(lambda,n,"lambda",1)), 
                                "assigning variable log_lik");
                }
            } else {

                current_statement_begin__ = 68;
                for (int n = 1; n <= n_row; ++n) {

                    current_statement_begin__ = 69;
                    stan::model::assign(log_lik, 
                                stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                                neg_binomial_2_log(get_base1(yint,n,"yint",1),get_base1(lambda,n,"lambda",1),get_base1(nb2_phi,1,"nb2_phi",1)), 
                                "assigning variable log_lik");
                }
            }

            // validate generated quantities
            current_statement_begin__ = 62;

            // write generated quantities
            for (int k_0__ = 0; k_0__ < n_row; ++k_0__) {
            vars__.push_back(log_lik[k_0__]);
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_bycatch";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (time_varying * (n_year - 1)); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "est_time_dev" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= time_varying; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_rw" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= family; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nb2_phi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pred" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= (time_varying * n_year); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "time_dev" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
        }


        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (time_varying * (n_year - 1)); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "est_time_dev" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= time_varying; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_rw" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= family; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nb2_phi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lambda" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pred" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_0__ = 1; k_0__ <= (time_varying * n_year); ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "time_dev" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
        }


        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= n_row; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }

}; // model

}

typedef model_bycatch_namespace::model_bycatch stan_model;


#endif