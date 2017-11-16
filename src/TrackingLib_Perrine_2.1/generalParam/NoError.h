/*******************************************************************************************************************
    This file defines the number error
    Institut Curie
    UMR - 144
    by Victor Racine
    28 october 2003
*******************************************************************************************************************/

#ifndef __NO_ERROR__
#define __NO_ERROR__


#define		RET_OK					0

//MIA
#define		SIZE_ERR				-1
#define		ALLOC_ERR				-2
#define		DIM_ERR					-3
#define		BOUND_ERR				-4
#define		WRITING_FILE_ERR		-5
#define		OPENING_FILE_ERR		-6
#define		READING_FILE_ERR		-7
#define		BAD_TYPE_READING_ERR	-8
#define		OPERATION_IMPOSSIBLE	-9
#define		DLL_ERROR				-10
#define		INIT_ERR				-11


//Segmentation Static
#define S_ALLOC_ERR					-201
#define S_OPEN_FILE_ERR				-202
#define S_WRTIE_ERR					-203
#define S_READ_ERR					-204

//morphoFilter Static
#define ERR_MORPHO_PARAM			-301
#define ERR_MORPHO_DIM				-302

//FFTStatic
#define ERR_FFT_ALLOC				-401
#define ERR_FFT_DIM					-402
#define ERR_FFT_SIZE				-403
#define ERR_FFT						-404

//ReadWriteSTK
#define ERR_READ_STK				-501
#define ERR_WRITE_STK				-502
#define ERR_BAD_FORMAT				-503
#define ERR_CLOSE_STK				-504
#define ERR_GRAPH_DIM				-505
#define ERR_GRAPH_SIZE				-506
#define ERR_TAG_NOT_FOUND			-507

//MIA MT
#define ERR_PARAM_MT				-601
#define ERR_GLOBAL_MT				-602

//RealtionMM_Tracking
#define	ERR_OPENNING_FILE			-701
#define ERR_EOF_UNEXPECTED			-702
#define ERR_READ_NB_NODES			-703
#define ERR_READ_NODES				-704
#define ERR_ALLOCATING				-705
#define ERR_FILE_FORMAT				-706


//Tracking
#define	T_ERR_OPENNING_FILE			-801
#define T_ERR_EOF_UNEXPECTED		-802
#define T_ERR_READ_NB_NODES			-803
#define T_ERR_READ_INPUT_PARAM		-804
#define T_ERR_NBDIM					-805
#define T_ERR_ALLOC					-806
#define T_ERR_CORRUPT_FILE			-807
#define T_ERR_BAD_MODEL				-808
#define T_ERR_OOP					-809
#define T_ERR_INVALIDE_DATA			-811

//Control of allocation
#define ERR_CONTROL_ALLOC_1			-901
#define ERR_CONTROL_ALLOC_2			-902
#define ERR_CONTROL_ALLOC_3			-903


//Misc Tools
#define ERR_DIM						-1001
#define ERR_UNKNOWN					-1002 //in execMIA

#endif //__NO_ERROR__
