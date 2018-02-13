#ifndef __MPC_H__
#define __MPC_H__

#include "libs.h"
#include "utillites.hpp"

class FG_eval
{
public:
	typedef CppAD::vector<CppAD::AD<double>> ADvector;
	void operator()(ADvector &fg, const ADvector &vars)
	{
		fg[0] = 0.0;
		for (uint32_t i = 0; i < data.N; ++i)
		{
			fg[0] += data.costs[0] * CppAD::pow(vars[data.starts[4] + i], 2);
			fg[0] += data.costs[1] * CppAD::pow(vars[data.starts[5] + i], 2);
			fg[0] += data.costs[2] * CppAD::pow(vars[data.starts[3] + i] - data.v_ref, 2);
		}
		for (uint32_t i = 0; i <  data.N - 1; ++i)
		{
			fg[0] += data.costs[3] * CppAD::pow(vars[data.starts[6] + i], 2);
			fg[0] += data.costs[4] * CppAD::pow(vars[data.starts[7] + i], 2);
			fg[0] += data.costs[7] * CppAD::pow(vars[data.starts[6] + i] * vars[data.starts[3] + i], 2);
		}
		for (uint32_t i = 0; i <  data.N - 2; ++i)
		{
			fg[0] += data.costs[5] * CppAD::pow(vars[data.starts[6] + i + 1] - vars[data.starts[6] + i], 2);
			fg[0] += data.costs[6] * CppAD::pow(vars[data.starts[7] + i + 1] - vars[data.starts[7] + i], 2);
		}
		for (uint32_t i = 0; i < data.starts.size() - 2; i++)
			fg[1 + data.starts[i]] = vars[data.starts[i]];
	
		for (uint32_t t = 1; t <  data.N; ++t)
		{
			CppAD::AD<double> x1		  = vars[data.starts[0] + t];
			CppAD::AD<double> y1		  = vars[data.starts[1] + t];
			CppAD::AD<double> psi1		  = vars[data.starts[2] + t];
			CppAD::AD<double> v1		  = vars[data.starts[3] + t];
			CppAD::AD<double> cte1		  = vars[data.starts[4] + t];
			CppAD::AD<double> epsi1		  = vars[data.starts[5] + t];

			CppAD::AD<double> x0		  = vars[data.starts[0] + t - 1];
			CppAD::AD<double> y0		  = vars[data.starts[1] + t - 1];
			CppAD::AD<double> psi0		  = vars[data.starts[2] + t - 1];
			CppAD::AD<double> v0		  = vars[data.starts[3] + t - 1];
			CppAD::AD<double> cte0		  = vars[data.starts[4] + t - 1];
			CppAD::AD<double> epsi0		  = vars[data.starts[5] + t - 1];
			CppAD::AD<double> delta	  	  = vars[data.starts[6] + t - 1];
			CppAD::AD<double> a		      = vars[data.starts[7] + t - 1];
			
			CppAD::AD<double> f0		  = coeff[0] + coeff[1] * x0 + coeff[2] * CppAD::pow(x0, 2) + coeff[3] * CppAD::pow(x0, 3);
			CppAD::AD<double> psides0	  = CppAD::atan(coeff[1] + 2 * coeff[2] * x0 + 3 * coeff[3] * CppAD::pow(x0, 2));
		
			fg[1 + data.starts[0] + t]    = x1    - (x0 + v0 * CppAD::cos(psi0) * data.dt);
			fg[1 + data.starts[1] + t] 	  = y1    - (y0 + v0 * CppAD::sin(psi0) * data.dt);	
			fg[1 + data.starts[2] + t] 	  = psi1  - (psi0 - v0/data.Lf * delta * data.dt);
			fg[1 + data.starts[3] + t] 	  = v1    - (v0 + a * data.dt);
			fg[1 + data.starts[4] + t] 	  = cte1  - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * data.dt));
			fg[1 + data.starts[5] + t] 	  = epsi1 - ((psi0 - psides0) - v0/data.Lf * delta * data.dt);
		}
	}

	void set_constants(const Data &_data)
	{
		data = _data;
	}
	void set_coeffs(const Eigen::VectorXd &_coeffs)
	{
		coeff = _coeffs;
	}

private:
	Data		    data;
	Eigen::VectorXd coeff;
};

struct MPC
{
	typedef CppAD::vector<double> Dvector;

	MPC();

	std::vector<double> solve(const Eigen::VectorXd &state, Eigen::VectorXd &coeffs);

	uint32_t get_polyfitorder() const;
	uint32_t get_state_size()   const;
	double get_Lf_constant()    const;
	double get_dt_constant() 	const;

private:
	void resize_vectors();
	void set_vectors();
	void set_starts();
	void set_solver_options();
	std::vector<double> set_result();

	CppAD::ipopt::solve_result<Dvector> solution;

	std::string					 options;

	FG_eval						 fg_eval;

	Dvector						 vars;
	Dvector						 vars_lowerbound;
	Dvector						 vars_upperbound;
	Dvector						 constraints_lowerbound;
	Dvector						 constraints_upperbound;
	
	bool						 solver_status;

	double						 dt 				= 0.1;
	double						 str_vars_upbound 	= 0.436332;
	double 						 str_vars_lowbound  = -0.436332;
	double						 thr_vars_upbound 	= 1.0;
	double 						 thr_vars_lowbound  = -1.0;
	double						 constr_lowbound 	= 0.0;
	double 						 constr_upbound 	= 0.0;
	double						 Lf 				= 2.67;
	double						 v_ref 				= 90;

	uint32_t					 N 					= 10;
	uint32_t					 n_state 			= 6;
	uint32_t					 n_control 			= 2;
	uint32_t					 n_vars;
	uint32_t					 n_constraints;
	uint32_t					 poly_fit_order 	= 3;
	
	std::vector<uint32_t>		 starts;
	std::vector<double>			 costs;
};
#endif