# Walkthrough: https://ropenscilabs.github.io/drake-manual/intro.html
# Download the code: drake_example("main") # nolint

# Load your packages and supporting functions into your session.
# If you use supporting scripts like the ones below,
# you will need to supply them yourself. Examples:
                                        # https://github.com/wlandau/drake-examples/tree/master/main/R
source("code/packages.R")  # Load your packages, e.g. library(drake).
source("code/config.R")
source("code/functions.R") # Define your custom code as a bunch of functions.
source("code/plan.R")      # Create your drake plan.
                                        # Call make() to run your work.
                                        # Your targets will be stored in a hidden .drake/ cache,
                                        # and you can read them back into memory with loadd() and read().




make(plan, parallelism = data_config$parallelism,
     jobs = data_config$jobs,
     prework = quote(future::plan(future::multisession)))
## make(plan)
