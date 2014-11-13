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
    {                   \
        LOG_ERROR(msg); \
        throw(-1);      \
    }

#define ASSERT_ALLWAYS(expression, msg)                   \
    {                                                     \
        if ((expression) == true)                         \
            ABORT(#expression " does not hold! " << msg); \
    }

#endif /* defined(__Fermion_Simulator__error_handler__) */
