#include "MPC.h"

MPC::MPC() 
{
	set_solver_options();
	
	n_vars		  = N *  n_state + (N - 1) * n_control;
	n_constraints = n_state*N;

	resize_vectors();
	set_starts();
	set_vectors();

	fg_eval.set_constants(Data(starts, costs, dt, v_ref, Lf, N));
}

std::vector<double> MPC::solve(const Eigen::VectorXd &state, Eigen::VectorXd &coeffs)
{
	for (uint32_t i = 0; i < starts.size() - 2; ++i)
	{
		vars[starts[i]]					  = state[i];
		constraints_lowerbound[starts[i]] = state[i];
		constraints_upperbound[starts[i]] = state[i];
	}
	
	fg_eval.set_coeffs(coeffs);

	CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
										  constraints_upperbound, fg_eval, solution);

	solver_status &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

	return set_result();
}

uint32_t MPC::get_polyfitorder() const
{
	return poly_fit_order;
}
uint32_t MPC::get_state_size()   const
{
	return n_state;
}
double MPC::get_Lf_constant() const
{
	return Lf;
}
double MPC::get_dt_constant() 	const
{
	return dt;
}


std::vector<double> MPC::set_result()
{
	std::vector<double> result;
	result.push_back(solution.x[starts[6]]);
	result.push_back(solution.x[starts[7]]);

	for (uint32_t i = 0; i < N - 1; ++i)
	{
		result.push_back(solution.x[starts[0] + i + 1]);
		result.push_back(solution.x[starts[1] + i + 1]);
	}
	return result;
}
void MPC::resize_vectors()
{
	vars.resize(n_vars);
	starts.resize(n_state + n_control);

	vars_lowerbound.resize(n_vars);
	vars_upperbound.resize(n_vars);

	constraints_lowerbound.resize(n_constraints);
	constraints_upperbound.resize(n_constraints);
	
	costs.resize(8);
}
void MPC::set_vectors()
{
	for (uint32_t i = 0; i < starts[starts.size() - 1]; ++i)
	{
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}
	for (uint32_t i = starts[starts.size() - 1]; i < starts.back(); ++i)
	{
		vars_lowerbound[i] = str_vars_lowbound;
		vars_upperbound[i] = str_vars_upbound;
	}
	for (uint32_t i = starts.back(); i < n_vars; ++i)
	{
		vars_lowerbound[i] = thr_vars_lowbound;
		vars_upperbound[i] = thr_vars_upbound;
	}
	for (uint32_t i = 0; i < n_constraints; ++i)
	{
		constraints_lowerbound[i] = constr_lowbound;
		constraints_upperbound[i] = constr_upbound;
	}
	
	costs[0] = 3500;
	costs[1] = 3500;
	costs[2] = 1;
	costs[3] = 5;
	costs[4] = 5;
	costs[5] = 200;
	costs[6] = 10;
	costs[7] = 10;
	
}
void MPC::set_starts()
{
	for (uint32_t i = 0; i < starts.size(); ++i)
	{
		if (i == 0)
			starts[0] = 0;
		else if (i == starts.size() - 1)
			starts[i] = starts[i - 1] + N - 1;
		else
			starts[i] = starts[i - 1] + N;
	}
}
void MPC::set_solver_options()
{
	options += "Integer print_level  0\n";
	options += "Sparse  true        forward\n";
	options += "Sparse  true        reverse\n";
	options += "Numeric max_cpu_time          0.5\n";
}