// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_BASE_H_
#define SRC_EOS_BASE_H_

#include "pretty_ostream.h"
#include "grid.h"

#include <string>
#include <vector>
#include <gsl/gsl_spline.h>

class EOS_base {
 private:
    int whichEOS;
    int number_of_tables;
    double eps_max;
    bool flag_muB;
    bool flag_muS;
    bool flag_muC;

 public:
    pretty_ostream music_message;
    std::vector<double> nb_bounds;
    std::vector<double> e_bounds;
    std::vector<double> fugl_bounds;
    std::vector<double> fugs_bounds;

    std::vector<double> nb_spacing;
    std::vector<double> e_spacing;
    std::vector<double> fugl_spacing;
    std::vector<double> fugs_spacing;

    std::vector<int> nb_length;
    std::vector<int> e_length;
    std::vector<int> fugl_length;
    std::vector<int> fugs_length;

    double ***pressure_tb;
    double ***temperature_tb;
    double ***mu_B_tb;
    double ***mu_S_tb;
    double ***mu_C_tb;
    double ***cs2_tb;
    double ***energy_tb;
    double ***entropy_tb;

    gsl_interp_accel *acc;
    gsl_spline *pressure_spline;
    gsl_spline *temperature_spline;
    gsl_spline *cs2_spline;
    gsl_spline *entropy_spline;

    EOS_base() = default;
    virtual ~EOS_base();

    std::string get_hydro_env_path() const;

    void set_number_of_tables(int ntables) {number_of_tables = ntables;}
    int  get_number_of_tables() const {return(number_of_tables);}
    void resize_table_info_arrays();

    void set_EOS_id(int eos_id) {whichEOS = eos_id;}
    int  get_EOS_id() const {return(whichEOS);}

    void set_flag_muB(bool flag_muB_in) {flag_muB = flag_muB_in;}
    bool get_flag_muB() const {return(flag_muB);}
    void set_flag_muS(bool flag_muS_in) {flag_muS = flag_muS_in;}
    bool get_flag_muS() const {return(flag_muS);}
    void set_flag_muC(bool flag_muC_in) {flag_muC = flag_muC_in;}
    bool get_flag_muC() const {return(flag_muC);}

    // returns maximum local energy density of the EoS table
    // in the unit of [1/fm^4]
    void   set_eps_max(double eps_max_in) {eps_max = eps_max_in;}
    double get_eps_max() const {return(eps_max);}

    double interpolate1D(double e, int table_idx, double ***table) const;
    double interpolate2D(double e, double rhob, int table_idx, double ***table) const;
    double interpolate3D(double e, double fugl, double fugs, double ***table) const;

    int    get_table_idx(double e) const;
    double get_entropy  (double epsilon, double rhob, double proper_tau) const;

    double calculate_velocity_of_sound_sq(double e, double rhob, double proper_tau) const;
    double get_dpOverde3(double e, double rhob, double proper_tau) const;
    double get_dpOverdrhob2(double e, double rhob, double proper_tau) const;
    double get_s2e_finite_rhob(double s, double rhob, double proper_tau) const;
    double get_T2e_finite_rhob(const double T, const double rhob, double proper_tau) const;

    virtual void   initialize_eos () {}
    virtual void   initialize_eos (int eos_id_in) {}
    virtual double get_cs2        (double e, double rhob, double proper_tau) const;
    virtual double p_rho_func     (double e, double rhob, double proper_tau) const {return(0.0);}
    virtual double p_e_func       (double e, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_temperature(double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_muB        (double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_muS        (double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_muC        (double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_rhoS       (double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_rhoC       (double epsilon, double rhob, double proper_tau) const {return(0.4*rhob);}
    virtual double get_pressure   (double epsilon, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_s2e        (double s, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_T2e        (double T, double rhob, double proper_tau) const {return(0.0);}
    virtual double get_light_fugacity (double proper_tau) const {return(0.0);}
    virtual double get_strange_fugacity (double proper_tau) const {return(0.0);}
    virtual void   check_eos      () const {}
    virtual double get_T_ratio    (Cell_small &cell) const {return(1.0);}
    
    void check_eos_with_finite_muB() const;
    void check_eos_no_muB() const;
};


#endif  // SRC_EOS_BASE_H_
