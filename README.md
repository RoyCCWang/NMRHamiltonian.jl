# NMRHamiltonian.jl
Given chemical shifts and J-coupling values, produce simulated amplitude and frequencies, and partitions of them.

See `/examples/simulate.jl` for an example.

# Install
From a Julia REPL or script, add the custom registry before adding the package:

```
using Pkg

Pkg.Registry.add(url = "https://github.com/RoyCCWang/RWPublicJuliaRegistry")

Pkt.add("NMRHamiltonian)
```

# Citation
Our work is undergoing peer review.

# License
This project is licensed under the GPL V3.0 license; see the LICENSE file for details. Individual source files may contain the following tag instead of the full license text:
```
SPDX-License-Identifier: GPL-3.0-only
```

Using SPDX enables machine processing of license information based on the SPDX License Identifiers and makes it easier for developers to see at a glance which license they are dealing with.
