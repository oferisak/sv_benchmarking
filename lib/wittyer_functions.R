# Functions to run wittyer
get_wittyer_command<-function(test_vcf,reference_vcf,output_folder,config_file=NA,target_bed=NA,em='SimpleCounting'){
  wittyer_command <-glue('docker run --workdir $(pwd) -v $(pwd):$(pwd) wittyer -i {test_vcf} -t {reference_vcf} --em "{em}" -o {output_folder}')
  if (!is.na(config_file)){
    wittyer_command<-paste0(wittyer_command,
                            glue(' --configFile {config_file}'),
                            collapse=' ')
  }
  if (!is.na(target_bed)){
    wittyer_command<-paste0(wittyer_command,
                            glue(' --includeBed {target_bed}'),
                            collapse=' ')
  }
  return(wittyer_command)
}

run_wittyer<-function(wittyer_command){
  command_output<-system(wittyer_command,intern = T)
  return(command_output)
}

parse_wittyer_ouput<-function(wittyer_output_folder){
  wittyer_stats_file<-glue('{wittyer_output_folder}/Wittyer.Stats.json')
  # Give the input file name to the function.
  result <- fromJSON(file = wittyer_stats_file)
  # Convert JSON file to a data frame.
  overall_stats<-data.frame(variant_type='All',bin='All',result$PerSampleStats[[1]]$OverallStats%>%list_bind_all())
  detailed_stats<-result$PerSampleStats[[1]]$DetailedStats
  variant_types<-c()
  per_bin_stats_df<-NULL
  for (vt in 1:length(detailed_stats)){variant_types<-c(variant_types,detailed_stats[[vt]]$VariantType)}
  names(detailed_stats)<-variant_types
  for (variant_type in variant_types){
    variant_overall_stats<-detailed_stats[[variant_type]]$OverallStats%>%list_bind_all()
    overall_stats<-overall_stats%>%
      bind_rows(data.frame(variant_type=variant_type,bin='All',variant_overall_stats))
    bin_names<-c()
    variant_all_bin_stats<-detailed_stats[[variant_type]]$PerBinStats
    for (bn in 1:length(variant_all_bin_stats)){bin_names<-c(bin_names,variant_all_bin_stats[[bn]]$Bin)}
    names(variant_all_bin_stats)<-bin_names
    
    for (bn in bin_names){
      bin_stats<-variant_all_bin_stats[[bn]]$Stats%>%list_bind_all()
      overall_stats<-overall_stats%>%bind_rows(data.frame(variant_type=variant_type,
                                                        bin=as.character(bn),
                                                        bin_stats))
    }
  }
  overall_stats<-overall_stats%>%mutate(across(where(is.list),function(x){unlist(x)}),
                                        bin=factor(bin,levels=c('[1, 50)',
                                                                '[1, 1000)',
                                                                '[50, 100)',
                                                                '[100, 200)',
                                                                '[200, 500)',
                                                                '[500, 1000)',
                                                                '[1000, 5000)',
                                                                '[5000, 10000)',
                                                                '10000+',
                                                                '[10000, 20000)',
                                                                '[20000, 50000)',
                                                                '[50, 1000)',
                                                                '1000+',
                                                                '50000+',
                                                                'All',
                                                                'NA'
                                                                )))
  
  overall_stats_list<-list('base'=overall_stats%>%filter(StatsType=='Base'),
                           'event'=overall_stats%>%filter(StatsType=='Event'))
  # write output tables
  write.table(overall_stats_list$event,
              file=glue('{wittyer_output_folder}/event_stats.tsv'),
              sep='\t',
              row.names = F)
  write.table(overall_stats_list$base,
              file=glue('{wittyer_output_folder}/base_stats.tsv'),
              sep='\t',
              row.names = F)
  return(overall_stats_list)
}

plot_wittyer_precision_recall<-function(wittyer_output_parsed,output_folder,base_or_event='event'){
  wittyer_output_events<-wittyer_output_parsed[[base_or_event]]
  g<-wittyer_output_events%>%select(variant_type,bin,Recall,Precision,QueryTotalCount)%>%
    filter(variant_type%in%c('Deletion','Insertion'),bin!='NA')%>%
    pivot_longer(-c(variant_type,bin,QueryTotalCount))%>%
    ggplot(aes(x=bin,y=as.numeric(value),size=QueryTotalCount,label=round(as.numeric(value),2)))+
    geom_point()+
    geom_text(hjust=-0.5, vjust=0,size=3)+
    facet_grid(variant_type~name,scales='free')+
    labs(x='Event size',y='',size='Number of events')+
    theme_minimal()+
    theme(legend.position = 'bottom',plot.background = element_rect(color='white'))
  print(g)
  ggsave(g,filename = glue('{output_folder}/per_bin_{base_or_event}_precision_recall_stats.jpg'),
         device = 'jpg',dpi = 150,height = 6,width = 12)
}
