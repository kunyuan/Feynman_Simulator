//
//  error_handler.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__error_handler__
#define __Fermion_Simulator__error_handler__

#include "logger.h"

#define ABORT(msg)      \
    do {                \
        LOG_ERROR(msg); \
        throw(-1);      \
    } while (0)

#define ASSERT_ALLWAYS(expression, msg)                   \
    do {                                                  \
        if ((expression) == false)                        \
            ABORT(#expression " does not hold! " << msg); \
    } while (0)

#endif /* defined(__Fermion_Simulator__error_handler__) */
