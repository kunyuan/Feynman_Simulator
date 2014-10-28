//
//  scopeguard.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/27/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_scopeguard_h
#define Feynman_Simulator_scopeguard_h

#define SCOPEGUARD_LINENAME_CAT(name, line) name##line
#define SCOPEGUARD_LINENAME(name, line) SCOPEGUARD_LINENAME_CAT(name, line)
#define ON_SCOPE_EXIT(callback) ScopeGuard SCOPEGUARD_LINENAME(EXIT, __LINE__)(callback)

class ScopeGuard {
  public:
    explicit ScopeGuard(std::function<void()> onExitScope)
        : onExitScope_(onExitScope), dismissed_(false)
    {
    }

    ~ScopeGuard()
    {
        if (!dismissed_) {
            onExitScope_();
        }
    }

    void Dismiss()
    {
        dismissed_ = true;
    }

  private:
    std::function<void()> onExitScope_;
    bool dismissed_;

  private: // noncopyable
    ScopeGuard(ScopeGuard const &);
    ScopeGuard &operator=(ScopeGuard const &);
};

#endif
