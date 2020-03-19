globus endpoint search "INTERVAL RNAseq humgen pipeline" # get identifier for endpoint
sangerRNA=7df826aa-20d9-11ea-9707-021304b0cca7 # set human-readable name
globus endpoint show $sangerRNA

globus ls $sangerRNA # list directory contents

