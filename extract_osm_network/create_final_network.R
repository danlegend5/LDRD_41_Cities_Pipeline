# Imports
library(data.table)
library(rgdal)
library(bit64)
library(igraph)
library(Imap)

# Clean cache
rm(list=ls())
options(stringsAsFactors=F)


#####################################################################################################
# Check the strong connectivity of the graph (i.e. can one travel from any node on the graph to any
# other node via the edges)
#####################################################################################################

# Read in the edge list and the node list
country_name = "australia"
city_name = "melbourne"
edgelist <- fread(paste0("/data/dmb20/Traffic.Data/Loop.Detector.Reformatted.Data.41.Cities/Results/Underlying.Network/", country_name, "/", city_name, "/", city_name, ".edg.csv"), fill=T, integer64="character")[edge_id!=""]
nodelist <- fread(paste0("/data/dmb20/Traffic.Data/Loop.Detector.Reformatted.Data.41.Cities/Results/Underlying.Network/", country_name, "/", city_name, "/", city_name, ".nod.csv"), integer64="character")

# Retain certain columns
edgelist <- edgelist[edge_to!="",.(edge_from,edge_to,edge_type,edge_id)]

# Join of the edge list and the node list
edgelist <- edgelist[nodelist,on=.(edge_from=node_id)]
edgelist <- edgelist[nodelist,on=.(edge_to=node_id)]
edgelist <- edgelist[!(is.na(node_x)|is.na(node_y)|is.na(i.node_x)|is.na(i.node_y))]

# Compute the length of each edge
edgelist[,L:=gdist(node_x,node_y,i.node_x,i.node_y,units = "m")]
edgelist <- edgelist[,.(edge_from,edge_to,L,edge_type,edge_id)]

# Create a graph
netgraph <- graph_from_data_frame(edgelist)

# Compute edge betweenness centrality
betw <- estimate_edge_betweenness(netgraph,cutoff=45)

# Copy and put in a table
edgelist_s <- edgelist[,.(edge_from,edge_to,edge_id,edge_type)]
betw_d <- data.table(betw=betw, edgelist_s, igraphid=as.numeric(E(netgraph)))

# Select residential roads with betweenness rank which is higher than 0.3
setorder(betw_d,-betw)
betw_d[edge_type=="highway.living_street", edge_type:="highway.residential"]
betw_d[!is.infinite(betw),rank:=1:.N,by=edge_type]
betw_d[,rel_rank:=rank/max(rank,na.rm = T),by=edge_type]
del_edges <- betw_d[(edge_type=="highway.residential" & rel_rank<0.3)]
del_edges[grep("-",edge_id),neg_edge:=gsub("-","",edge_id)]
del_edges[!(grep("-",edge_id)),neg_edge:=paste0("-",edge_id)]
lookup <- unique(betw_d[,.(e_id=edge_id, i_id=igraphid)])
del_edges[lookup,neg_igraphid:=i_id,on=.(neg_edge=e_id)]
ilookup <- data.table(edges=E(netgraph),numericid=1:length(E(netgraph)))
udel <- unique(c(del_edges$igraphid,del_edges$neg_igraphid))

# Check strong connectedness
netgraph_mod <- delete_edges(netgraph,ilookup[numericid %in% udel,edges])
ccstong <- components(netgraph_mod, mode = c("strong"))
sel <- data.table(csize=ccstong$csize)
sel[,id:=1:.N]
subgr <- sel[csize==max(csize)]
subs <- data.table(ccstong$membership,names(ccstong$membership))
subs <- subs[V1==subgr$id]
edges <- edgelist[,.(edge_from,edge_to,edge_id)]
edges <- edges[(edge_from %in% subs$V2) & (edge_to %in% subs$V2)]
edges_valid <- data.table(edges$edge_id)

# Write out the edges we want to keep
write.table(edges_valid, file = paste0("/data/dmb20/Traffic.Data/Loop.Detector.Reformatted.Data.41.Cities/Results/Underlying.Network/", country_name, "/", city_name, "/", "edges_valid_strongcon.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)


#####################################################################################################
# RUN NETCONVERT IN SUMO and remove the edges not in the file
#####################################################################################################
#
# shell(paste0("cd ..\\Sumo\\",folder," && netconvert -s net/osm.net.xml --keep-edges.input-file edges_valid_strongcon.txt -o net/",city,"_osm.net.xml --proj.plain-geo true --geometry.remove true --plain-output-prefix ",city," -v true"))
#
#####################################################################################################
# MODIFY INTERSECTIONS
#####################################################################################################
# CHANGE rightbeforeleft to allWay stops allway_stop (USA) in node file!
# CHANGE tls to tls turn on red traffic_light_right_on_red (USA) in node file!
# # shell(paste0("cd ..\\Sumo\\RIO14A && netconvert -s RIO_osm.net.xml  --proj.plain-geo true --plain-output-prefix RIO -v true"))
# shell(paste0("python C:\\Users\\ambuehll\\Documents\\sumo-msvc12extrax64-git\\tools\\xml\\xml2csv.py ","C:/Users/ambuehll/Sumo/",folder,"/RIO.edg.xml -s \";\""))
# shell(paste0("python C:\\Users\\ambuehll\\Documents\\sumo-msvc12extrax64-git\\tools\\xml\\xml2csv.py ","C:/Users/ambuehll/Sumo/",folder,"/RIO.nod.xml --xsd C:/Users/ambuehll/Sumo/",folder,"/xsd/nodes_file.xsd"))
#
# node <- fread(paste0("C:/Users/ambuehll/Sumo/",folder,"/",city,".nod.csv"))
# edge <- fread(paste0("C:/Users/ambuehll/Sumo/",folder,"/",city,".edg.csv"))
#
# nodes <- node[,.(node_id,node_type)]
#
# # -----------------------
# # STEP 1:
#   # ZIPPER FINDER
#   # -----------------------
#
#   tt <- unique(edge[,.(edge_from,edge_to,edge_id,edge_priority,edge_type,edge_numLanes)])
#   unique(tt[,.(edge_type,edge_priority)])
#
#   # change links to real links for prio rules in zipper
#
#   tt[edge_type=="highway.trunk_link",edge_priority:=13]
#   tt[edge_type=="highway.trunk_link",edge_type:="highway.trunk"]
#
#   leaving <- tt[,.(N_leave=length(unique(edge_id))),edge_from]
#   tt <- tt[,.(list(unique(edge_priority)),N=length(unique(sub("-","",sub("#.*","",edge_id)))),list(unique(edge_type)), list(unique(edge_numLanes))),by=edge_to]
#
#   tt[,Le:=sapply(tt$V1, length)]
#   tt[,unumLan:=sapply(tt$V4, length)]
#
#   tt <- tt[leaving,on=c(edge_to="edge_from")]
#
#   interest <- tt[N==2 & Le==1 & N_leave==1 & unumLan==1]
#
#   interest <- interest[V3 %in% c("highway.motorway","highway.motorway_link","highway.primary","highway.secondary","highway.trunk")]
#   interesting <- nodes[node_id %in% interest$edge_to]
#
#   thatx <- interesting[!(node_type %in% c("traffic_light_right_on_red","traffic_light","allway_stop","dead_end"))]
#   node[node_id %in% thatx$node_id,node_type:="zipper"]
#
#   # -----------------------
#   # STEP 2
#   # CALCULATE ADDITIONAL TRAFFIC SIGNALS
#   # -----------------------
#
#   # STEP 2a)
#   # intersections with a high number of lanes and high priorities coming together should be signal controlled!
#   # only for OSM network with low quality!
#   if(city=="RIO"){
#
#   tt <- edge[,.(edge_from,edge_to,edge_id,edge_priority,edge_type,edge_numLanes)]
#   tt[edge_type=="highway.tertiary",edge_priority:=7L]
#   leaving <- tt[,.(N_leave=length(unique(edge_id))),edge_from]
#
#   tt <- unique(tt[,.(edge_to,edge_id,edge_type,edge_priority,edge_numLanes)])
#
#   tt <- tt[,.(list((edge_priority*edge_numLanes)),N=length(unique(sub("-","",sub("#.*","",edge_id)))),list(edge_type)),by=edge_to]
#
#   tt[,Le:=sapply(tt$V1, length)]
#   tt <- tt[leaving,on=c(edge_to="edge_from")]
#
#   interest <- tt[N>1]
#   interest <- nodes[interest,on=c(node_id="edge_to")]
#   motor <- interest[,.(paste0(unlist(V3))),by=.(node_id)]
#   # exclude motorways
#   motor <- motor[,motorway:=V1 %like% "motor"]
#
#   motor <- motor[,.(SM=sum(motorway)),by=node_id]
#   interest[motor,motor:=SM,on="node_id"]
#
#   interest[,sumPrio:=sum(unlist(V1)),by=node_id]
#   interesting <- interest[sumPrio >= 55 & node_type!= "traffic_light" & node_type!= "zipper" & node_type!= "dead_end" & motor==0 & N_leave>1]
#   node[node_type=="priority" & node_id %in% unique(interesting$node_id),node_type:="traffic_light"]
#   }
#
#   # ------
#   # STEP 2b)
#   # complicated intersections should be traffic controlled!
#   # creates a buffer around each intersection and evaluates how many overlap
#   # ------
#   devtools::install_github("valentinitnelav/geobuffer",force = T)
#   library(geobuffer)
#   library(tidyr)
#
#   roundabout_nodes <- edge[roundabout_nodes!=""]
#   roundabout_nodes <- separate_rows(roundabout_nodes, roundabout_nodes)
#   roundabout_nodes <- unique(roundabout_nodes$roundabout_nodes)
#
#   tt <- node[!(node_type %in% c("dead_end","traffic_light","zipper")),.(node_id,node_x,node_y,node_type)]
#   tt <- tt[!(node_id %in% roundabout_nodes)]
#   edg_nodes <- edge[!(edge_type %in% c("highway.residential","highway.living_street", "highway.unclassified")),.(edge_to,edge_from,edge_type)]
#   edg_nodes <- edg_nodes[edge_to!=""]
#   highwaynodes <- edg_nodes[edge_type %in% c("highway.trunk","highway.trunk_link","highway.motorway","highway.motorway_link")]
#   highwaynodes <- unique(c(highwaynodes$edge_from, highwaynodes$edge_to))
#   tt <- tt[!(node_id %in% highwaynodes)]
#   pts <- tt[,.(node_x,node_y)]
#
#   pts_buf <- geobuffer_pts(xy = pts, dist_m = 10)
#   tz <- data.table(over(pts_buf,pts_buf,returnList=T))
#   tt[,closest:=tz$V1]
#   tt[,Le:=sapply(closest, length)]
#   complicated_intersections <- tt[Le>3,.(closest)]
#   complicated_intersections <- unique(separate_rows(complicated_intersections,sep = ","))
#   complicated <- tt[complicated_intersections$closest,]
#
#   node[node_id %in% unique(complicated$node_id),node_type:="traffic_light"]
#
# # node[node_type=="traffic_light" ,node_type:="traffic_light_right_on_red"]
# # node[node_type=="right_before_left" ,node_type:="allway_stop"]
# node <- node[,.(node_id,node_x,node_y,node_type,node_tl)]
#
# #     <node id="1016407632" x="105855.80" y="92027.22" type="traffic_light_right_on_red" tl="1016407632"/>
# cat("<nodes>\n", file=paste0("C:/Users/ambuehll/Sumo/",folder,"/",city",N.xml"),append=F)
#
# for(i in 1:nrow(node)){
#
#   # print(i)
#   if(node[i]$node_tl==""){
#     tt <- paste0("<node id=\"",node[i]$node_id ,"\" x= \"",node[i]$node_x,"\" y= \"",node[i]$node_y,"\" type= \"",node[i]$node_type,"\" />\n")
#   }else{
#     tt <- paste0("<node id=\"",node[i]$node_id ,"\" x= \"",node[i]$node_x,"\" y= \"",node[i]$node_y,"\" type= \"",node[i]$node_type,"\" tl= \"",node[i]$node_tl,"\"  />\n")
#   }
#
#   cat(tt, file="C:/Users/ambuehll/Sumo/RIO14A/RION.xml",append=T)
#
# }
#
# cat("</nodes>", file="C:/Users/ambuehll/Sumo/RIO14A/RION.xml",append=T)
#
# # ----
# shell(paste0("cd ..\\Sumo\\RIO14A && netconvert -n RION.xml -e RIO.edg.xml -t RIO.typ.xml -o RIO.net.xml --tls.join true --tls.uncontrolled-within true --proj.plain-geo true --proj.utm true --no-internal-links true --no-turnarounds.tls true --no-turnarounds.geometry true --plain-output-prefix RIOT --geometry.remove false -v true"))
# # ###################################################################################################################
# # Load new edge and node files
# # Create TAZ and Demand files
# # ###################################################################################################################
#
# shell(paste0("python C:\\Users\\ambuehll\\Documents\\sumo-msvc12extrax64-git\\tools\\xml\\xml2csv.py ","C:/Users/ambuehll/Sumo/RIO14A/RIO.edg.xml"))
# shell(paste0("python C:\\Users\\ambuehll\\Documents\\sumo-msvc12extrax64-git\\tools\\xml\\xml2csv.py ","C:/Users/ambuehll/Sumo/RIO14A/RIO.nod.xml"))
#
# # delete last cols in excel!
# edgelist <- fread("C:/Users/ambuehll/Sumo/RIO14A/RIO.edg.csv",fill=T,integer64 = "character")[edge_id!=""]
# nodelist <- fread("C:/Users/ambuehll/Sumo/RIO14A/RIO.nod.csv",integer64 = "character")
#
# shp <- readOGR("C:/Users/ambuehll/Sumo/RIO14A/tract/RIO.shp")
# raster::extent(shp)
# CRS(proj4string(shp))
# sumonodes <- SpatialPointsDataFrame(nodelist[,.(node_x,node_y)], nodelist, coords.nrs = numeric(0),
#                                     proj4string = CRS(proj4string(shp)), match.ID=F)
#
# library(maptools)
# grants.districts <- over(sumonodes,shp)  # Get district data
# sumonodes <- spCbind(sumonodes, grants.districts)
# sumonodes <- data.table(sumonodes@data)
#
# edgelist <- edgelist[!is.null(edge_from),.(edge_from,edge_id,edge_to)]
# edgelist <- edgelist[sumonodes,on=.(edge_from=node_id)]
#
# tazasco <- edgelist[,.(taz=as.numeric(Zona),gid=edge_id)]
#
# # create TAZ
#
# tazasco <- tazasco[!is.na(gid)]
# tazasco <- tazasco[!is.na(taz)]
#
# # test flows
# # ODS[tazasco,vorh:=1,on=.(o=taz)]
# # ODS[,sum(flow),by=vorh]
#
# setorder(tazasco,taz,gid)
# tazasco <- tazasco[!is.na(taz)]
# tazasco[,iter:=.GRP,by=taz]
#
# # <taz id="1" edges="c(1616, 1630, 13827, 21612, 30598, 30820)"/>
# cat("<tazs>\n", file="C:/Users/ambuehll/Sumo/RIO14A/taz.xml",append=F)
#
# for(i in 1:max(tazasco$iter)){
#
#   print(i)
#   su <- paste0(tazasco[iter==i]$gid,collapse = " ")
#   tt <- paste0("<taz id=\"",unique(tazasco[iter==i]$taz) ,"\" edges= \"",su,"\"/>\n")
#
#   cat(tt, file="C:/Users/ambuehll/Sumo/RIO14A/taz.xml",append=T)
#
# }
#
# cat("</tazs>", file="C:/Users/ambuehll/Sumo/RIO14A/taz.xml",append=T)
#
# # points(sumonodes,add=T,col="blue")
# # create sumo
#
# taz <- unique(tazasco$taz)
# # taz <- fread("C:/Users/ambuehll/Sumo/RIO_osm_small/taz.csv")$taz_id
#
# lifi <- list.files("C:/Users/ambuehll/Sumo/RIO14A/vehicle_trips_rect/",
#                    pattern = "",full.names = T)
#
# # -------------------------------------
# # ODS[,flow:=all_flow]
# # ODS[,o:=origin_tract]
# # ODS[,d:=dest_tract]
# tt <- list()
# for(i in lifi){
#   ODS <- fread(i)
#   ODS <- ODS[(o %in% taz) & (d %in% taz)]
#   timeint <- substr(x = i, 64,65)
#   tt[[i]] <- data.table(ODS[,sum(flow)])
#   print(ODS[,sum(flow)])
#
# cat("$OR;D2
# * From-Time  To-Time\n",
#       file=paste0("C:/Users/ambuehll/Sumo/RIO14A/OD/OD",as.numeric(timeint),".fma"))
#
#   write(paste0(paste0(as.numeric(timeint),".","00")," ",paste0(as.numeric(timeint)+1,".","00")),
#         file=paste0("C:/Users/ambuehll/Sumo/RIO14A/OD/OD",as.numeric(timeint),".fma"),append=T)
#
# cat("* Factor
# 1.00
# *\n", file=paste0("C:/Users/ambuehll/Sumo/RIO14A/OD/OD",as.numeric(timeint),".fma"),append=T)
#
#   write.table(ODS[,.(o,d,flow)], paste0("C:/Users/ambuehll/Sumo/RIO14A/OD/OD",as.numeric(timeint),".fma"), sep = " ", dec = ".",
#               row.names = F,append=T, col.names = F)
#
# }
#
# # ##############################################################
# # CREATE DEMAND
# # ##############################################################
#
# shell(paste0("cd ..\\Sumo\\RIO14A && OD2trips -n taz.xml -s 1 -d \"OD\\OD3.fma,OD\\OD4.fma,OD\\OD5.fma,OD\\OD6.fma,OD\\OD7.fma,OD\\OD8.fma,OD\\OD9.fma,OD\\OD10.fma,OD\\OD11.fma,OD\\OD12.fma\" -o trips24h.trips.xml --departlane best --departpos random_free --different-source-sink true"))
