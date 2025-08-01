#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void odin_model_emanators_initmod_desolve(void *);
extern void odin_model_emanators_output_dde(void *);
extern void odin_model_emanators_rhs_dde(void *);
extern void odin_model_emanators_rhs_desolve(void *);
extern void odin_model_hrp2_initmod_desolve(void *);
extern void odin_model_hrp2_output_dde(void *);
extern void odin_model_hrp2_rhs_dde(void *);
extern void odin_model_hrp2_rhs_desolve(void *);
extern void odin_model_initmod_desolve(void *);
extern void odin_model_IVM_SMChet_initmod_desolve(void *);
extern void odin_model_IVM_SMChet_output_dde(void *);
extern void odin_model_IVM_SMChet_rhs_dde(void *);
extern void odin_model_IVM_SMChet_rhs_desolve(void *);
extern void odin_model_mass_effect_initmod_desolve(void *);
extern void odin_model_mass_effect_output_dde(void *);
extern void odin_model_mass_effect_pp_initmod_desolve(void *);
extern void odin_model_mass_effect_pp_output_dde(void *);
extern void odin_model_mass_effect_pp_rhs_dde(void *);
extern void odin_model_mass_effect_pp_rhs_desolve(void *);
extern void odin_model_mass_effect_rhs_dde(void *);
extern void odin_model_mass_effect_rhs_desolve(void *);
extern void odin_model_output_dde(void *);
extern void odin_model_rhs_dde(void *);
extern void odin_model_rhs_desolve(void *);
extern void odin_model_TBV_initmod_desolve(void *);
extern void odin_model_TBV_output_dde(void *);
extern void odin_model_TBV_rhs_dde(void *);
extern void odin_model_TBV_rhs_desolve(void *);

/* .Call calls */
extern SEXP odin_model_contents(SEXP);
extern SEXP odin_model_create(SEXP);
extern SEXP odin_model_emanators_contents(SEXP);
extern SEXP odin_model_emanators_create(SEXP);
extern SEXP odin_model_emanators_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_emanators_metadata(SEXP);
extern SEXP odin_model_emanators_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_emanators_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_emanators_set_user(SEXP, SEXP);
extern SEXP odin_model_hrp2_contents(SEXP);
extern SEXP odin_model_hrp2_create(SEXP);
extern SEXP odin_model_hrp2_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_hrp2_metadata(SEXP);
extern SEXP odin_model_hrp2_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_hrp2_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_hrp2_set_user(SEXP, SEXP);
extern SEXP odin_model_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_IVM_SMChet_contents(SEXP);
extern SEXP odin_model_IVM_SMChet_create(SEXP);
extern SEXP odin_model_IVM_SMChet_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_IVM_SMChet_metadata(SEXP);
extern SEXP odin_model_IVM_SMChet_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_IVM_SMChet_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_IVM_SMChet_set_user(SEXP, SEXP);
extern SEXP odin_model_mass_effect_contents(SEXP);
extern SEXP odin_model_mass_effect_create(SEXP);
extern SEXP odin_model_mass_effect_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_mass_effect_metadata(SEXP);
extern SEXP odin_model_mass_effect_pp_contents(SEXP);
extern SEXP odin_model_mass_effect_pp_create(SEXP);
extern SEXP odin_model_mass_effect_pp_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_mass_effect_pp_metadata(SEXP);
extern SEXP odin_model_mass_effect_pp_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_mass_effect_pp_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_mass_effect_pp_set_user(SEXP, SEXP);
extern SEXP odin_model_mass_effect_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_mass_effect_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_mass_effect_set_user(SEXP, SEXP);
extern SEXP odin_model_metadata(SEXP);
extern SEXP odin_model_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_set_user(SEXP, SEXP);
extern SEXP odin_model_TBV_contents(SEXP);
extern SEXP odin_model_TBV_create(SEXP);
extern SEXP odin_model_TBV_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_TBV_metadata(SEXP);
extern SEXP odin_model_TBV_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_TBV_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_TBV_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"odin_model_emanators_initmod_desolve",      (DL_FUNC) &odin_model_emanators_initmod_desolve,      1},
    {"odin_model_emanators_output_dde",           (DL_FUNC) &odin_model_emanators_output_dde,           1},
    {"odin_model_emanators_rhs_dde",              (DL_FUNC) &odin_model_emanators_rhs_dde,              1},
    {"odin_model_emanators_rhs_desolve",          (DL_FUNC) &odin_model_emanators_rhs_desolve,          1},
    {"odin_model_hrp2_initmod_desolve",           (DL_FUNC) &odin_model_hrp2_initmod_desolve,           1},
    {"odin_model_hrp2_output_dde",                (DL_FUNC) &odin_model_hrp2_output_dde,                1},
    {"odin_model_hrp2_rhs_dde",                   (DL_FUNC) &odin_model_hrp2_rhs_dde,                   1},
    {"odin_model_hrp2_rhs_desolve",               (DL_FUNC) &odin_model_hrp2_rhs_desolve,               1},
    {"odin_model_initmod_desolve",                (DL_FUNC) &odin_model_initmod_desolve,                1},
    {"odin_model_IVM_SMChet_initmod_desolve",     (DL_FUNC) &odin_model_IVM_SMChet_initmod_desolve,     1},
    {"odin_model_IVM_SMChet_output_dde",          (DL_FUNC) &odin_model_IVM_SMChet_output_dde,          1},
    {"odin_model_IVM_SMChet_rhs_dde",             (DL_FUNC) &odin_model_IVM_SMChet_rhs_dde,             1},
    {"odin_model_IVM_SMChet_rhs_desolve",         (DL_FUNC) &odin_model_IVM_SMChet_rhs_desolve,         1},
    {"odin_model_mass_effect_initmod_desolve",    (DL_FUNC) &odin_model_mass_effect_initmod_desolve,    1},
    {"odin_model_mass_effect_output_dde",         (DL_FUNC) &odin_model_mass_effect_output_dde,         1},
    {"odin_model_mass_effect_pp_initmod_desolve", (DL_FUNC) &odin_model_mass_effect_pp_initmod_desolve, 1},
    {"odin_model_mass_effect_pp_output_dde",      (DL_FUNC) &odin_model_mass_effect_pp_output_dde,      1},
    {"odin_model_mass_effect_pp_rhs_dde",         (DL_FUNC) &odin_model_mass_effect_pp_rhs_dde,         1},
    {"odin_model_mass_effect_pp_rhs_desolve",     (DL_FUNC) &odin_model_mass_effect_pp_rhs_desolve,     1},
    {"odin_model_mass_effect_rhs_dde",            (DL_FUNC) &odin_model_mass_effect_rhs_dde,            1},
    {"odin_model_mass_effect_rhs_desolve",        (DL_FUNC) &odin_model_mass_effect_rhs_desolve,        1},
    {"odin_model_output_dde",                     (DL_FUNC) &odin_model_output_dde,                     1},
    {"odin_model_rhs_dde",                        (DL_FUNC) &odin_model_rhs_dde,                        1},
    {"odin_model_rhs_desolve",                    (DL_FUNC) &odin_model_rhs_desolve,                    1},
    {"odin_model_TBV_initmod_desolve",            (DL_FUNC) &odin_model_TBV_initmod_desolve,            1},
    {"odin_model_TBV_output_dde",                 (DL_FUNC) &odin_model_TBV_output_dde,                 1},
    {"odin_model_TBV_rhs_dde",                    (DL_FUNC) &odin_model_TBV_rhs_dde,                    1},
    {"odin_model_TBV_rhs_desolve",                (DL_FUNC) &odin_model_TBV_rhs_desolve,                1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"odin_model_contents",                          (DL_FUNC) &odin_model_contents,                          1},
    {"odin_model_create",                            (DL_FUNC) &odin_model_create,                            1},
    {"odin_model_emanators_contents",                (DL_FUNC) &odin_model_emanators_contents,                1},
    {"odin_model_emanators_create",                  (DL_FUNC) &odin_model_emanators_create,                  1},
    {"odin_model_emanators_initial_conditions",      (DL_FUNC) &odin_model_emanators_initial_conditions,      2},
    {"odin_model_emanators_metadata",                (DL_FUNC) &odin_model_emanators_metadata,                1},
    {"odin_model_emanators_rhs_r",                   (DL_FUNC) &odin_model_emanators_rhs_r,                   3},
    {"odin_model_emanators_set_initial",             (DL_FUNC) &odin_model_emanators_set_initial,             4},
    {"odin_model_emanators_set_user",                (DL_FUNC) &odin_model_emanators_set_user,                2},
    {"odin_model_hrp2_contents",                     (DL_FUNC) &odin_model_hrp2_contents,                     1},
    {"odin_model_hrp2_create",                       (DL_FUNC) &odin_model_hrp2_create,                       1},
    {"odin_model_hrp2_initial_conditions",           (DL_FUNC) &odin_model_hrp2_initial_conditions,           2},
    {"odin_model_hrp2_metadata",                     (DL_FUNC) &odin_model_hrp2_metadata,                     1},
    {"odin_model_hrp2_rhs_r",                        (DL_FUNC) &odin_model_hrp2_rhs_r,                        3},
    {"odin_model_hrp2_set_initial",                  (DL_FUNC) &odin_model_hrp2_set_initial,                  4},
    {"odin_model_hrp2_set_user",                     (DL_FUNC) &odin_model_hrp2_set_user,                     2},
    {"odin_model_initial_conditions",                (DL_FUNC) &odin_model_initial_conditions,                2},
    {"odin_model_IVM_SMChet_contents",               (DL_FUNC) &odin_model_IVM_SMChet_contents,               1},
    {"odin_model_IVM_SMChet_create",                 (DL_FUNC) &odin_model_IVM_SMChet_create,                 1},
    {"odin_model_IVM_SMChet_initial_conditions",     (DL_FUNC) &odin_model_IVM_SMChet_initial_conditions,     2},
    {"odin_model_IVM_SMChet_metadata",               (DL_FUNC) &odin_model_IVM_SMChet_metadata,               1},
    {"odin_model_IVM_SMChet_rhs_r",                  (DL_FUNC) &odin_model_IVM_SMChet_rhs_r,                  3},
    {"odin_model_IVM_SMChet_set_initial",            (DL_FUNC) &odin_model_IVM_SMChet_set_initial,            4},
    {"odin_model_IVM_SMChet_set_user",               (DL_FUNC) &odin_model_IVM_SMChet_set_user,               2},
    {"odin_model_mass_effect_contents",              (DL_FUNC) &odin_model_mass_effect_contents,              1},
    {"odin_model_mass_effect_create",                (DL_FUNC) &odin_model_mass_effect_create,                1},
    {"odin_model_mass_effect_initial_conditions",    (DL_FUNC) &odin_model_mass_effect_initial_conditions,    2},
    {"odin_model_mass_effect_metadata",              (DL_FUNC) &odin_model_mass_effect_metadata,              1},
    {"odin_model_mass_effect_pp_contents",           (DL_FUNC) &odin_model_mass_effect_pp_contents,           1},
    {"odin_model_mass_effect_pp_create",             (DL_FUNC) &odin_model_mass_effect_pp_create,             1},
    {"odin_model_mass_effect_pp_initial_conditions", (DL_FUNC) &odin_model_mass_effect_pp_initial_conditions, 2},
    {"odin_model_mass_effect_pp_metadata",           (DL_FUNC) &odin_model_mass_effect_pp_metadata,           1},
    {"odin_model_mass_effect_pp_rhs_r",              (DL_FUNC) &odin_model_mass_effect_pp_rhs_r,              3},
    {"odin_model_mass_effect_pp_set_initial",        (DL_FUNC) &odin_model_mass_effect_pp_set_initial,        4},
    {"odin_model_mass_effect_pp_set_user",           (DL_FUNC) &odin_model_mass_effect_pp_set_user,           2},
    {"odin_model_mass_effect_rhs_r",                 (DL_FUNC) &odin_model_mass_effect_rhs_r,                 3},
    {"odin_model_mass_effect_set_initial",           (DL_FUNC) &odin_model_mass_effect_set_initial,           4},
    {"odin_model_mass_effect_set_user",              (DL_FUNC) &odin_model_mass_effect_set_user,              2},
    {"odin_model_metadata",                          (DL_FUNC) &odin_model_metadata,                          1},
    {"odin_model_rhs_r",                             (DL_FUNC) &odin_model_rhs_r,                             3},
    {"odin_model_set_initial",                       (DL_FUNC) &odin_model_set_initial,                       4},
    {"odin_model_set_user",                          (DL_FUNC) &odin_model_set_user,                          2},
    {"odin_model_TBV_contents",                      (DL_FUNC) &odin_model_TBV_contents,                      1},
    {"odin_model_TBV_create",                        (DL_FUNC) &odin_model_TBV_create,                        1},
    {"odin_model_TBV_initial_conditions",            (DL_FUNC) &odin_model_TBV_initial_conditions,            2},
    {"odin_model_TBV_metadata",                      (DL_FUNC) &odin_model_TBV_metadata,                      1},
    {"odin_model_TBV_rhs_r",                         (DL_FUNC) &odin_model_TBV_rhs_r,                         3},
    {"odin_model_TBV_set_initial",                   (DL_FUNC) &odin_model_TBV_set_initial,                   4},
    {"odin_model_TBV_set_user",                      (DL_FUNC) &odin_model_TBV_set_user,                      2},
    {NULL, NULL, 0}
};

void R_init_ICDMM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
