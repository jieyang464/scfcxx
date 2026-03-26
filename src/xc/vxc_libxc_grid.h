#pragma once

#include "xc/vxc_evaluator.h"

namespace xc {

enum BuiltinFunctionalId {
    None = 0,
    LsdaExchange = 1,
};

}  // namespace xc

VxcFunctor make_libxc_vxc_functor();
