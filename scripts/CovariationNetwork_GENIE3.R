library(magrittr)
library(tidyr)
library(data.table)

# PMS treatment network

# this is the expression matrix including all modalities where each row was labelled ModalityName_ID
exp.mat = read.csv('analysis/multiOmics_expression_matrix.csv',row.names = 1)

# table of curated pathways and entities
curated.selections = read.csv('data/multiOmics_curated_pathways.csv')
curated.selections = tidyr::pivot_longer(curated.selections,cols = 3:ncol(curated.selections))
curated.selections$value = trimws(curated.selections$value)

curated.selections = curated.selections %>% 
  separate_rows(value, sep=",")
curated.selections$displayed_name = paste0(curated.selections$name,'_',curated.selections$value)
curated.selections = curated.selections[curated.selections$value!='',]
colnames(curated.selections)[1:2]=c('Hub','HubTwo')

# only keep the proteins/metabolites that are related to the pathways of interest
exp.mat = exp.mat[rownames(exp.mat)%in%curated.selections$displayed_name,]

set.seed(2023)
res <- GENIE3::GENIE3(as.matrix(exp.mat), targets = rownames(exp.mat),
                      nCores = 5)
saveRDS(res,'data/multiOmics_GENIE3Output.rds')
#res = readRDS('data/multiOmics_GENIE3Output.rds')
links <- bulkAnalyseR::get_link_list_rename(res,20)

# function from bulkAnalyseR
get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) %>%
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

plot_GRN_multiomics <- function(weightMat, plotConnections){
  
  # get the normal edges and nodes from GENIE3
  edges <- get_link_list_rename(weightMat, plotConnections)
  # this bit removes duplicate edges where just the start and end are reversed
  edges$id = paste0(pmin(edges$from,edges$to),pmax(edges$from,edges$to))
  edges = edges[order(-edges$value),]
  edges = edges[!duplicated(edges$id),]
  edges = edges[,1:3]
  # adding edge colours
  edges$color = ifelse(sub("\\_.*", "", edges$from)=='IntraLip','#B32F71',
                       ifelse(sub("\\_.*", "", edges$from)=='Secretomics','#117733',
                              ifelse(sub("\\_.*", "", edges$from)=='IntraProt','#0DBB47',
                                     ifelse(sub("\\_.*", "", edges$from)=='ExtraMet','#332288',
                                            ifelse(sub("\\_.*", "", edges$from)=='IntraMet','#6E5DC3','darkgreen')))))
  
  # add nodes for each end of the edges and colour them by modality
  nodes <- tibble::tibble(
    id = c(edges$to, edges$from),
    label = '',
    group = sub("\\_.*", "", .data$id),
    color = ifelse(sub("\\_.*", "", .data$id)=='IntraLip','#B32F71',
                   ifelse(sub("\\_.*", "", .data$id)=='Secretomics','#117733',
                          ifelse(sub("\\_.*", "", .data$id)=='IntraProt','#0DBB47',
                                 ifelse(sub("\\_.*", "", .data$id)=='ExtraMet','#332288',
                                        ifelse(sub("\\_.*", "", .data$id)=='IntraMet','#6E5DC3','darkgreen'))))),
    value=1,
    font.color='white',
    font.size=30,
    font.bold=F
  ) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  
  # pathway edges and nodes
  # add edges between all entities and all hub three nodes if they exist, otherwise the hub two node
  pathway.edges = data.frame('from'=curated.selections$HubTwo,
                             'to'=curated.selections$displayed_name,
                             'value'=min(edges$value),
                             'color'='#BFE4F7')
  pathway.nodes.hub2 = data.frame('id'=unique(curated.selections$HubTwo),
                                  'label'=unique(curated.selections$HubTwo),
                                  'group'='Hub2',
                                  'color'='#BFE4F7',
                                  'value'=1,
                                  'font.color'='#36454F',
                                  'font.size'=800,
                                  'font.bold'=F)
  pathway.nodes.hub = data.frame('id'=unique(curated.selections$Hub),
                                 'label'=unique(curated.selections$Hub),
                                 'group'='Hub',
                                 'color'='#BFE4F7',
                                 'value'=1,
                                 'font.color'='black',
                                 'font.size'=1500,
                                 'font.bold'=T)
  pathway.edges.hub2.hub = data.frame('from'=curated.selections$HubTwo,
                                      'to'=curated.selections$Hub,
                                      'value'=3*max(edges$value),
                                      'color'='black')
  # combine all the pathways edges and nodes
  pathway.edges = rbind(pathway.edges,pathway.edges.hub2.hub)
  pathway.edges = unique(pathway.edges)
  pathway.nodes = rbind(pathway.nodes.hub,pathway.nodes.hub2)
  pathway.edges = pathway.edges[pathway.edges$from!=pathway.edges$to,]
  # combine the original edges and nodes with the pathway edges and nodes
  nodes = rbind(nodes,pathway.nodes) 
  edges = rbind(edges,pathway.edges)
  edges = edges[edges$from != edges$to,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  edges = edges[edges$to %in% nodes$id,]
  edges = edges[edges$from %in% nodes$id,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  # add line breaks for long strings
  nodes$label = gsub('_','\n',nodes$label)
  nodes$label = gsub('/','/\n',nodes$label)
  nodes$label = ifelse(nodes$group!='Hub',stringr::str_wrap(nodes$label,width = 10),
                       stringr::str_wrap(nodes$label,width = 100))
  nodes = nodes[nodes$id!='',]
  
  # built the network using the edges and nodes you've defined
  # the physics part depends how busy your network is
  network <- visNetwork::visNetwork(nodes, edges,width=3000,height=1500) %>%
    visNetwork::visGroups(groupname = "IntraLip", color = '#B32F71', shape = "triangle") %>%
    visNetwork::visGroups(groupname = "Secretomics", color = '#117733', shape = "dot") %>%
    visNetwork::visGroups(groupname = "Hub2", color = '#BFE4F7', shape = "box") %>%
    visNetwork::visGroups(groupname = "Hub", color = '#BFE4F7', shape = "box") %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           forceAtlas2Based = list(gravitationalConstant = -4000,springConstant=0.01,centralGravity=0.005)) %>%
    visNetwork::visLegend()%>%
    visNetwork::visLayout(randomSeed = 2023)
  return(list('nodes'=nodes,'edges'=edges,'network'=network))
}

# output the network and save to html
plot_GRN_multiomics(res,1000)
visNetwork::visSave(plot_GRN_multiomics(res,1000)$network, 'analysis/GENIE3_multiOmics_1000.html', selfcontained = TRUE, background = "white")
