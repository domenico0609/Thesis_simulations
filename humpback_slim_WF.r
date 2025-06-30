library(slimr)

##Slim script for the simulation of the nuclear and mitochondrial DNA
slim_script(
  slim_block(initialize(), 
             {
               ##initializing a nucleotide based non-WF simulation
               initializeSLiMModelType("WF");	
               initializeSLiMOptions(nucleotideBased=T);
               
               ##reading a Fasta file as Ancestral string for the individuals
               defineConstant("L", initializeAncestralNucleotides("C:/Users/Domenico/Documents/humpy/Ancestral_Fasta_test.fasta"));
               
               ##Sexual reproduction call
               initializeSex("A");
               
               ##Length of the mitochondrial section
               defineConstant("M", 21);
               
               ##number of populations
               defineConstant("P", 4);
               
               ##Number of total individuals in a single population
               defineConstant("Nmax", 200);
               
               ##Number of founding individuals 
               defineConstant("N", 10);
               
               ##Number of generations
               defineConstant("Gen", 10000);
               
               for (i in 1:10) initializeMutationTypeNuc(i, 0.5, "f", 0.0);
               initializeMutationType("m11", 0.5, "f", 0.0);
               
               m1.color = "purple";
               m2.color = "pink";
               m3.color = "blue";
               m4.color = "red";
               m5.color = "yellow";
               m6.color = "green";
               m7.color = "blue";
               m8.color = "red";
               m9.color = "yellow";
               m10.color = "green";
               m11.color = "grey";
               
               initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(asFloat(0.0)));	##neutral nuclear DNA region
               initializeGenomicElementType("g2", m2, 1.0, mmJukesCantor(asFloat(0.0)));	##neutral mitochondrial DNA region 
               initializeGenomicElement(g1, 0, L-M-1);
               initializeGenomicElement(g2, L-M, L-1);
               initializeRecombinationRate(c(1e-8, 0), c(L-M-1, L-1));	
               
               g1.color = "blue";
               g2.color = "pink";	
               
             }
  ),
  
  slim_block(mateChoice(),
             {	
               matePop = sample(sim.subpopulations, 1);
               while(matePop == individual.subpopulation) matePop = sample(sim.subpopulations, 1);
               mate= matePop.sampleIndividuals(1, sex = "M");
               return(mate);
  }
  ),
  
  slim_block(1, early(),
             {
    for (i in 1:P) sim.addSubpop(i, N);
    
    p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1N.vcf", m3);
    p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1M.vcf", m7);
    p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2N.vcf", m4);	
    p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2M.vcf", m8);
    p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3N.vcf", m5);	
    p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3M.vcf", m9);
    p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4N.vcf", m6);	
    p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4M.vcf", m10);	
    
    log = community.createLogFile("sim_log.txt", logInterval=Gen);
    log.addCycle();
    log.addCustomColumn("FST(p1&p2)", "calcFST(p1.genomes,p2.genomes);");
    log.addCustomColumn("FST(p1&p3)", "calcFST(p1.genomes,p3.genomes);");	
    log.addCustomColumn("FST(p1&p4)", "calcFST(p1.genomes,p4.genomes);");
    log.addCustomColumn("FST(p2&p3)", "calcFST(p2.genomes,p3.genomes);");	
    log.addCustomColumn("FST(p2&p4)", "calcFST(p2.genomes,p4.genomes);");
    log.addCustomColumn("FST(p3&p4)", "calcFST(p3.genomes,p4.genomes);");
    
    gen1 = sim.subpopulations.individuals;
    gen1.genome1.addNewMutation(m11, 0.0, L-1);
    
    }
  ),
  
  slim_block(10, early(),
    {
        for (subpop in sim.subpopulations) subpop.setSubpopulationSize(Nmax);     
    }
  ),
  
  slim_block(modifyChild(), {
    
    if (!child.genome1.containsMarkerMutation(m11, L-1)){
      FMp1 = parent1.genome1.mutationsOfType(m7);
      FMp2 = parent1.genome1.mutationsOfType(m8);
      FMp3 = parent1.genome1.mutationsOfType(m9);
      FMp4 = parent1.genome1.mutationsOfType(m10);
      FMMark = parent1.genome1.mutationsOfType(m11);
      FMMut = parent1.genome1.mutationsOfType(m2);	
      
      child.genome1.addMutations(c(FMp1, FMp2, FMp3,FMp4, FMMark, FMMut));
    }
    
    CMp1 = child.genome2.mutationsOfType(m7);
    CMp2 = child.genome2.mutationsOfType(m8);
    CMp3 = child.genome2.mutationsOfType(m9);
    CMp4 = child.genome2.mutationsOfType(m10)
    CMMark = child.genome2.mutationsOfType(m11);
    CMMut = child.genome2.mutationsOfType(m2);
    
    child.genome2.removeMutations(c(CMp1, CMp2, CMp3, CMp4, CMMark, CMMut));
    
    return(T);
  }
  ),
  
  slim_block(1000, late(), { 
    
    ## output nucleotide data with slimr verb.
    r_output_nucleotides();
    r_output_nucleotides(both_genomes = TRUE);
    sim.outputFixedMutations(); } )) -> humpy_script

results <- slim_run(humpy_script)
res_data <- slim_results_to_data(results)

dev.new()
image(ape::as.DNAbin(res_data$data[[2]]))