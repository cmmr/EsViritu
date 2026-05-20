# Attribution for vendored dataui code

This directory contains a vendored subset of the **dataui** R package,
used to render genome coverage sparklines in EsViritu HTML reports.

**Source repository**: https://github.com/timelyportfolio/dataui  
**Upstream version**: 0.0.1 (master branch)  
**License**: MIT

## MIT License (dataui)

Copyright (c) 2020 Kent Russell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## What is vendored

- `R/dataui_sparklines.R`: sparkline-relevant R functions from `R/duisparkline.R`,
  `R/components.R`, `R/reactable_helpers.R`, `R/dependency.R`, and `R/utils.R`
- `inst/htmlwidgets/dataui.js`: compiled htmlwidget JS bundle
- `inst/www/dataui.standalone.js`: compiled standalone JS bundle

The JavaScript assets (`inst/`) are also covered by the MIT license and additionally
incorporate the **data-ui** library by Chris Williams
(https://williaster.github.io/data-ui), also MIT licensed.
