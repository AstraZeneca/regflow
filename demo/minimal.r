library( regflow )

import.database.omnipath("./database_omnipath.rda", TRUE)


interacs <- get.interactions( "omnipath", 
    "./database_omnipath.rda", "./table_interacs.csv",
    "./cache_interacs.rda", FALSE )

reg.graph <- build.regulation.graph( interacs, NULL, "./cache_reg_graph.rda", 
    FALSE )

reg.matrix <- calculate.regulation.matrices( reg.graph, 
    "./cache_reg_matrix.rda", FALSE )

write.regulation.flow( reg.matrix$flow, "MAPK3", FALSE, "." )

plot_regulation.flow( reg.matrix$flow, "PTPN11", "MAPK3", FALSE, ".", plot.seed = 12345 )

write.pathway( reg.matrix$flow, "PTPN11", "MAPK3", "." )

plot_pathway( reg.graph, reg.matrix$step, reg.matrix$flow, "PTPN11", "MAPK3", 
    "." )

