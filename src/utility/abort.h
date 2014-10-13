//
//  error_handler.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__error_handler__
#define __Fermion_Simulator__error_handler__

#define ABORT(msg)          \
    {                       \
        LOGGER(ERROR, msg); \
        throw(-1);          \
    }

#endif /* defined(__Fermion_Simulator__error_handler__) */
