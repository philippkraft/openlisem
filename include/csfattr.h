#ifndef CSF__ATTR_H
#define CSF__ATTR_H

#ifndef lint
#define RCS_ID_CSFLEGEN_H "$Header: c:/lisemroot/CVSROOT/lisemroot/libsrc/libcsf/include/csfattr.h,v 1.1.1.1 2003/04/08 09:34:25 jetten2 Exp $"
#endif

#ifdef __cplusplus
 extern "C" {
#endif

typedef enum CSF_ATTR_ID {
        ATTR_ID_LEGEND_V1=1,   /* version 1 legend */
        ATTR_ID_HISTORY=2,     /* history fields */
        ATTR_ID_COLOUR_PAL=3,  /* colour palette */
        ATTR_ID_GREY_PAL=4,    /* grey palette */
        ATTR_ID_DESCRIPTION=5, /* description */
        ATTR_ID_LEGEND_V2=6    /* version 2 legend */
} CSF_ATTR_ID;

#define CSF_LEGEND_ENTRY_SIZE  64
#define CSF_LEGEND_DESCR_SIZE  60

typedef struct CSF_LEGEND {
	INT4    nr;
	char    descr[60];
} CSF_LEGEND;

#ifdef __cplusplus
 }
#endif

#endif /*  CSF__ATTR_H */