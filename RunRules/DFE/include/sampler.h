/**\file */
#ifndef SLIC_DECLARATIONS_sampler_H
#define SLIC_DECLARATIONS_sampler_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define sampler_map_number (5)
#define sampler_N (128)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_samplerKernel The number of ticks for which kernel "samplerKernel" will run.
 * \param [in] inscalar_samplerKernel_Tftinv0 Input scalar parameter "samplerKernel.Tftinv0".
 * \param [in] inscalar_samplerKernel_Tftinv1 Input scalar parameter "samplerKernel.Tftinv1".
 * \param [in] inscalar_samplerKernel_Tftinv2 Input scalar parameter "samplerKernel.Tftinv2".
 * \param [in] inscalar_samplerKernel_Tftinv3 Input scalar parameter "samplerKernel.Tftinv3".
 * \param [in] inscalar_samplerKernel_Tftinv4 Input scalar parameter "samplerKernel.Tftinv4".
 * \param [in] inscalar_samplerKernel_gridsquare_inv Input scalar parameter "samplerKernel.gridsquare_inv".
 * \param [in] inscalar_samplerKernel_iteration_scalar Input scalar parameter "samplerKernel.iteration_scalar".
 * \param [in] inscalar_samplerKernel_tau0 Input scalar parameter "samplerKernel.tau0".
 * \param [in] inscalar_samplerKernel_tau1 Input scalar parameter "samplerKernel.tau1".
 * \param [in] inscalar_samplerKernel_tau2 Input scalar parameter "samplerKernel.tau2".
 * \param [in] inscalar_samplerKernel_tau3 Input scalar parameter "samplerKernel.tau3".
 * \param [in] inscalar_samplerKernel_tau4 Input scalar parameter "samplerKernel.tau4".
 * \param [in] instream_Nbar Stream "Nbar".
 * \param [in] instream_size_Nbar The size of the stream instream_Nbar in bytes.
 * \param [in] instream_Sft Stream "Sft".
 * \param [in] instream_size_Sft The size of the stream instream_Sft in bytes.
 * \param [in] instream_data Stream "data".
 * \param [in] instream_size_data The size of the stream instream_data in bytes.
 * \param [in] instream_normal_variates Stream "normal_variates".
 * \param [in] instream_size_normal_variates The size of the stream instream_normal_variates in bytes.
 * \param [in] instream_normal_variates2 Stream "normal_variates2".
 * \param [in] instream_size_normal_variates2 The size of the stream instream_normal_variates2 in bytes.
 * \param [in] instream_sandtDFE Stream "sandtDFE".
 * \param [in] instream_size_sandtDFE The size of the stream instream_sandtDFE in bytes.
 * \param [out] outstream_reconDFE Stream "reconDFE".
 * \param [in] outstream_size_reconDFE The size of the stream outstream_reconDFE in bytes.
 */
void sampler(
	uint64_t ticks_samplerKernel,
	double inscalar_samplerKernel_Tftinv0,
	double inscalar_samplerKernel_Tftinv1,
	double inscalar_samplerKernel_Tftinv2,
	double inscalar_samplerKernel_Tftinv3,
	double inscalar_samplerKernel_Tftinv4,
	double inscalar_samplerKernel_gridsquare_inv,
	uint64_t inscalar_samplerKernel_iteration_scalar,
	double inscalar_samplerKernel_tau0,
	double inscalar_samplerKernel_tau1,
	double inscalar_samplerKernel_tau2,
	double inscalar_samplerKernel_tau3,
	double inscalar_samplerKernel_tau4,
	const void *instream_Nbar,
	size_t instream_size_Nbar,
	const void *instream_Sft,
	size_t instream_size_Sft,
	const void *instream_data,
	size_t instream_size_data,
	const void *instream_normal_variates,
	size_t instream_size_normal_variates,
	const void *instream_normal_variates2,
	size_t instream_size_normal_variates2,
	const void *instream_sandtDFE,
	size_t instream_size_sandtDFE,
	void *outstream_reconDFE,
	size_t outstream_size_reconDFE);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] ticks_samplerKernel The number of ticks for which kernel "samplerKernel" will run.
 * \param [in] inscalar_samplerKernel_Tftinv0 Input scalar parameter "samplerKernel.Tftinv0".
 * \param [in] inscalar_samplerKernel_Tftinv1 Input scalar parameter "samplerKernel.Tftinv1".
 * \param [in] inscalar_samplerKernel_Tftinv2 Input scalar parameter "samplerKernel.Tftinv2".
 * \param [in] inscalar_samplerKernel_Tftinv3 Input scalar parameter "samplerKernel.Tftinv3".
 * \param [in] inscalar_samplerKernel_Tftinv4 Input scalar parameter "samplerKernel.Tftinv4".
 * \param [in] inscalar_samplerKernel_gridsquare_inv Input scalar parameter "samplerKernel.gridsquare_inv".
 * \param [in] inscalar_samplerKernel_iteration_scalar Input scalar parameter "samplerKernel.iteration_scalar".
 * \param [in] inscalar_samplerKernel_tau0 Input scalar parameter "samplerKernel.tau0".
 * \param [in] inscalar_samplerKernel_tau1 Input scalar parameter "samplerKernel.tau1".
 * \param [in] inscalar_samplerKernel_tau2 Input scalar parameter "samplerKernel.tau2".
 * \param [in] inscalar_samplerKernel_tau3 Input scalar parameter "samplerKernel.tau3".
 * \param [in] inscalar_samplerKernel_tau4 Input scalar parameter "samplerKernel.tau4".
 * \param [in] instream_Nbar Stream "Nbar".
 * \param [in] instream_size_Nbar The size of the stream instream_Nbar in bytes.
 * \param [in] instream_Sft Stream "Sft".
 * \param [in] instream_size_Sft The size of the stream instream_Sft in bytes.
 * \param [in] instream_data Stream "data".
 * \param [in] instream_size_data The size of the stream instream_data in bytes.
 * \param [in] instream_normal_variates Stream "normal_variates".
 * \param [in] instream_size_normal_variates The size of the stream instream_normal_variates in bytes.
 * \param [in] instream_normal_variates2 Stream "normal_variates2".
 * \param [in] instream_size_normal_variates2 The size of the stream instream_normal_variates2 in bytes.
 * \param [in] instream_sandtDFE Stream "sandtDFE".
 * \param [in] instream_size_sandtDFE The size of the stream instream_sandtDFE in bytes.
 * \param [out] outstream_reconDFE Stream "reconDFE".
 * \param [in] outstream_size_reconDFE The size of the stream outstream_reconDFE in bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *sampler_nonblock(
	uint64_t ticks_samplerKernel,
	double inscalar_samplerKernel_Tftinv0,
	double inscalar_samplerKernel_Tftinv1,
	double inscalar_samplerKernel_Tftinv2,
	double inscalar_samplerKernel_Tftinv3,
	double inscalar_samplerKernel_Tftinv4,
	double inscalar_samplerKernel_gridsquare_inv,
	uint64_t inscalar_samplerKernel_iteration_scalar,
	double inscalar_samplerKernel_tau0,
	double inscalar_samplerKernel_tau1,
	double inscalar_samplerKernel_tau2,
	double inscalar_samplerKernel_tau3,
	double inscalar_samplerKernel_tau4,
	const void *instream_Nbar,
	size_t instream_size_Nbar,
	const void *instream_Sft,
	size_t instream_size_Sft,
	const void *instream_data,
	size_t instream_size_data,
	const void *instream_normal_variates,
	size_t instream_size_normal_variates,
	const void *instream_normal_variates2,
	size_t instream_size_normal_variates2,
	const void *instream_sandtDFE,
	size_t instream_size_sandtDFE,
	void *outstream_reconDFE,
	size_t outstream_size_reconDFE);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_samplerKernel; /**<  [in] The number of ticks for which kernel "samplerKernel" will run. */
	double inscalar_samplerKernel_Tftinv0; /**<  [in] Input scalar parameter "samplerKernel.Tftinv0". */
	double inscalar_samplerKernel_Tftinv1; /**<  [in] Input scalar parameter "samplerKernel.Tftinv1". */
	double inscalar_samplerKernel_Tftinv2; /**<  [in] Input scalar parameter "samplerKernel.Tftinv2". */
	double inscalar_samplerKernel_Tftinv3; /**<  [in] Input scalar parameter "samplerKernel.Tftinv3". */
	double inscalar_samplerKernel_Tftinv4; /**<  [in] Input scalar parameter "samplerKernel.Tftinv4". */
	double inscalar_samplerKernel_gridsquare_inv; /**<  [in] Input scalar parameter "samplerKernel.gridsquare_inv". */
	uint64_t inscalar_samplerKernel_iteration_scalar; /**<  [in] Input scalar parameter "samplerKernel.iteration_scalar". */
	double inscalar_samplerKernel_tau0; /**<  [in] Input scalar parameter "samplerKernel.tau0". */
	double inscalar_samplerKernel_tau1; /**<  [in] Input scalar parameter "samplerKernel.tau1". */
	double inscalar_samplerKernel_tau2; /**<  [in] Input scalar parameter "samplerKernel.tau2". */
	double inscalar_samplerKernel_tau3; /**<  [in] Input scalar parameter "samplerKernel.tau3". */
	double inscalar_samplerKernel_tau4; /**<  [in] Input scalar parameter "samplerKernel.tau4". */
	const void *instream_Nbar; /**<  [in] Stream "Nbar". */
	size_t instream_size_Nbar; /**<  [in] The size of the stream instream_Nbar in bytes. */
	const void *instream_Sft; /**<  [in] Stream "Sft". */
	size_t instream_size_Sft; /**<  [in] The size of the stream instream_Sft in bytes. */
	const void *instream_data; /**<  [in] Stream "data". */
	size_t instream_size_data; /**<  [in] The size of the stream instream_data in bytes. */
	const void *instream_normal_variates; /**<  [in] Stream "normal_variates". */
	size_t instream_size_normal_variates; /**<  [in] The size of the stream instream_normal_variates in bytes. */
	const void *instream_normal_variates2; /**<  [in] Stream "normal_variates2". */
	size_t instream_size_normal_variates2; /**<  [in] The size of the stream instream_normal_variates2 in bytes. */
	const void *instream_sandtDFE; /**<  [in] Stream "sandtDFE". */
	size_t instream_size_sandtDFE; /**<  [in] The size of the stream instream_sandtDFE in bytes. */
	void *outstream_reconDFE; /**<  [out] Stream "reconDFE". */
	size_t outstream_size_reconDFE; /**<  [in] The size of the stream outstream_reconDFE in bytes. */
} sampler_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void sampler_run(
	max_engine_t *engine,
	sampler_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *sampler_run_nonblock(
	max_engine_t *engine,
	sampler_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void sampler_run_group(max_group_t *group, sampler_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *sampler_run_group_nonblock(max_group_t *group, sampler_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void sampler_run_array(max_engarray_t *engarray, sampler_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *sampler_run_array_nonblock(max_engarray_t *engarray, sampler_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* sampler_convert(max_file_t *maxfile, sampler_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* sampler_init(void);

/* Error handling functions */
int sampler_has_errors(void);
const char* sampler_get_errors(void);
void sampler_clear_errors(void);
/* Free statically allocated maxfile data */
void sampler_free(void);
/* These are dummy functions for hardware builds. */
int sampler_simulator_start(void);
int sampler_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_sampler_H */

