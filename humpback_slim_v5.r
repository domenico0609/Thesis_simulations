library(slimr)

##Slim script for the simulation of the nuclear and mitochondrial DNA
slim_script(
  slim_block(initialize(), 
    {
	  ##initializing a nucleotide based non-WF simulation
	  initializeSLiMModelType("nonWF");	
	  initializeSLiMOptions(nucleotideBased=T);
    
    ##reading a Fasta file as Ancestral string for the individuals
	  defineConstant("L", initializeAncestralNucleotides("C:/Users/Domenico/Documents/humpy/Ancestral_Fasta_test.fasta"));
		
    ##Sexual reproduction call
	  initializeSex("A");	
	
	  ##vector of carrying capacities for each subpopulation and a vector for the benefit of the specific feeding ground of the subpopulation
  	defineConstant("K", c(100.0, 100.0, 100.0));	
	  defineConstant("Feed", c(1.1, 1.1, 1.1));
	
	  ##number of populations
	  defineConstant("P", 3);
	
	  ##Number of individuals per subpopulation at the start of the simulation 
	  defineConstant("N", 10);	
	
	  ##Age of maturity and maximum age of the population are defined here (will be changed to a sex-based age of maturity later on. This would require additionial changed to the mortality age table
	  defineConstant("Age_mat", 10);
	  defineConstant("Age_max", 80);
	
	  ##to counteract the early boost in fitness from the carrying capacity formula a low-density value is defined so that at low density mortality is only due to age and sex.
	  defineConstant("Low_density_benefit", asFloat(1.0));
	
    #Sex based survivability
	  defineConstant("Sur_male", 0.94);
	  defineConstant("Sur_female", 0.90);
	
	  Mort = float(Age_max+1);
	  Mort[(Age_max-10):Age_max] = c(0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);
	  Mort[0:Age_mat] = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.025, 0.025, 0.025, 0.025);
	
	  defineConstant ("Mort_age", Mort);
	
	  for (i in 1:8) initializeMutationTypeNuc(i, 0.5, "f", 0.0);
	  initializeMutationType("m9", 0.5, "f", 0.0);

  	m1.color = "yellow";
	  m2.color = "red";
	  m3.color = "orange";
	  m4.color = "blue";
	  m5.color = "cyan";
	  m6.color = "green";
	  m7.color = "pink";
	  m8.color = "purple";
	  m9.color = "gray";
	
  	initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(asFloat(0.0)));	##neutral nuclear DNA region
	  initializeGenomicElementType("g2", m2, 1.0, mmJukesCantor(asFloat(0.0)));	##neutral mitochondrial DNA region 
	  initializeGenomicElement(g1, 0, asInteger(L/2-1));
	  initializeGenomicElement(g2, asInteger(L/2), L-1);
	  initializeRecombinationRate(c(1e-8, 0), c(asInteger(L/2-1), L-1));	
	
	  g1.color = "blue";
	  g2.color = "pink";	
	
    }
  ),
  
  slim_block(reproduction(), 
	
	  {	
	  if (individual.sex == "F" & individual.age >= Age_mat & individual.tag == 0) {
	
		  matePop = sample(sim.subpopulations, 1);
		  mate= matePop.sampleIndividuals(1, sex = "M", minAge = Age_mat);
		
		  ##Checks if the above sample actually contains a male
		  if (size(mate) > 0)
		  { 
		  child = subpop.addCrossed(individual, mate);
		  child.tag = 0;
		  }
		
		  individual.tag = rdunif(1, min=1, max=2);
	  }
	
	    ##if the individual has a tag higher than 0 it is lowered by 1 and once the tag reaches zero again the female can reproduce again
	    else if (individual.tag > 0) individual.tag = individual.tag - 1;
	  }
  ),

  slim_block(1, early(),
           {
	  for (i in 1:P) sim.addSubpop(i, N);
	
    p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1N.vcf", m3);
    p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1M.vcf", m6);
    p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2N.vcf", m4);	
    p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2M.vcf", m7);
    p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3N.vcf", m5);	
    p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3M.vcf", m8);
	
	  for (subpop in sim.subpopulations) 
	  {
	    ind = subpop.individuals;
	    ind.age = rdunif(N, min=0, max=Age_max-10);
	    ind.tag = 0;
	  }		
	
	  gen1 = sim.subpopulations.individuals;
	  gen1.genome1.addNewMutation(m9, 0.0, L-1);
	
    }
  ),

  slim_block(early(),
    {
	  o = 0;
		
		##life table mortality p1
	  for (subpop in sim.subpopulations)
	    {
		  if (subpop.individualCount > 0) {
			  inds = subpop.individuals;
			  ages = inds.age;
			  age_mortality = Mort_age[ages];
			  age_survival = 1 - age_mortality;
	
			  ##fitness is calculated for all individuals based on their age and feeding ground
			  inds.fitnessScaling = age_survival * Feed[o];
	
			  ##the fitnessScaling is altered by a survival factor based on the sex of the individual
			  females = 0;
			  males = 0;		
		
			  for (ind in inds)
				  if (ind.sex =="F"){
					  ind.fitnessScaling = ind.fitnessScaling * Sur_female;		
					  females = females + 1;
					  }	
				  else if (ind.sex =="M"){
					  ind.fitnessScaling = ind.fitnessScaling * Sur_male;	
				  	males = males + 1;
				  }
		
			  maleRatio = males / (males + females);
			  femaleRatio = 1 - maleRatio;
				
			  sexComp = (Sur_male * maleRatio + Sur_female * femaleRatio);
		
			  subpop.fitnessScaling = min(K[o] / min(subpop.individualCount *
			  mean(age_survival) * Feed [o] * sexComp) , Low_density_benefit);	
		
			
		  }
		  o = o + 1;
	  }
	
    }
  ),
  
  slim_block(modifyChild(), {

	  if (!child.genome1.containsMarkerMutation(m9, L-1)){
	  FMp1 = parent1.genome1.mutationsOfType(m6);
	  FMp2 = parent1.genome1.mutationsOfType(m7);
	  FMp3 = parent1.genome1.mutationsOfType(m8);
	  FMMark = parent1.genome1.mutationsOfType(m9);
	  FMMut = parent1.genome1.mutationsOfType(m2);	
	
	  child.genome1.addMutations(c(FMp1, FMp2, FMp3, FMMark, FMMut));
	  }
	
	  CMp1 = child.genome2.mutationsOfType(m6);
	  CMp2 = child.genome2.mutationsOfType(m7);
	  CMp3 = child.genome2.mutationsOfType(m8);
	  CMMark = child.genome2.mutationsOfType(m9);
	  CMMut = child.genome2.mutationsOfType(m2);
	
	  child.genome2.removeMutations(c(CMp1, CMp2, CMp3, CMMark, CMMut));
		
	  return(T);
    }
  ),

slim_block(2000, late(), { 
  
    ## output nucleotide data with slimr verb.
    r_output_nucleotides();
    r_output_nucleotides(both_genomes = TRUE);
    sim.outputFixedMutations(); } )) -> humpy_script

results <- slim_run(humpy_script)
res_data <- slim_results_to_data(results)


Nucleo1 <- Biostrings::DNAStringSet("A")
Mito1 <- Biostrings::DNAStringSet("A")
Nucleo2 <- Biostrings::DNAStringSet("A")
Mito2 <- Biostrings::DNAStringSet("A")

for (i in 1:length(res_data$data[[1]])) {
  Nucleo1[[i]] <- res_data$data[[1]][[i]][1:50]
}

for (i in 1:length(res_data$data[[2]])) {
  Nucleo2[[i]] <- res_data$data[[2]][[i]][1:50]
}

for (i in 1:length(res_data$data[[1]])) {
  Mito1[[i]] <- res_data$data[[1]][[i]][51:71]
}

for (i in 1:length(res_data$data[[2]])) {
  Mito2[[i]] <- res_data$data[[2]][[i]][51:71]
}

dev.new()
image(ape::as.DNAbin(res_data$data[[1]]))
dev.new()
image(ape::as.DNAbin(res_data$data[[2]]))
dev.new()
image(ape::as.DNAbin(Nucleo1))
dev.new()
image(ape::as.DNAbin(Mito1))
dev.new()
image(ape::as.DNAbin(Nucleo2))
dev.new()
image(ape::as.DNAbin(Mito2))
