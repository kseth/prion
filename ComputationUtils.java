import java.util.*;

class ComputationUtils {

	public static void main(String ... args) {

	}

	private static int samplePolymerLength(HashMap<Integer, Integer> polymerLengths, int polymers) {
		
		double r, maxR;
		maxR = 0;
		int chosenPolymer = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet()) {
			r = Math.random()*entry.value()*1.0/polymers;
			if(r > maxR) {
				maxR = r;
				chosenPolymer = entry.key();
			}
		}

		return chosenPolymer;
	}

	private static int samplePolymerBond(HashMap<Integer, Integer> polymerLengths, int polymerBonds) {
		double r, maxR;
		maxR = 0;
		int chosenPolymer = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet()) {
			r = Math.random()*(entry.value())*(entry.key()-1)*1.0/polymerBonds;
			if(r > maxR) {
				maxR = r;
				chosenPolymer = entry.key();
			}
		}
	}

	public static HashMap<Integer, Integer> prionGillespiePolymerLengths(double endTime, double lambda, double delta_m, double beta, double delta_p, double b, int polymerThreshold, HashMap<Integer, Integer> polymerLengths0, int monomers0) {
		
		// variables to keep track of state
		double currentTime = 0;
		HashMap<Integer, Integer> polymerLengths = (HashMap<Integer, Integer>)(polymerLengths0.clone());
		int polymers = 0;
		for(int i: polymerLengths.values())
			polymers+=i;
		int monomers = monomers0;
		int polymerSubunits = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
			polymerSubunits+=entry.key()*entry.value();
		int polymerBonds = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
			polymerBonds+=(entry.key()-1)*entry.value();


		// variables to keep track of probabilities
		double aLambda, aDeltaM, aBeta, aDeltaP, aB, a;
		double r1, r2, waitTime;
		int destroyedPolymer, addedToPolymer, breakingPolymer, breakPosition, largeFragment, smallFragment;

		while(currentTime < endTime) {
			r1 = Math.random();
			r2 = Math.random();

			aLambda = lambda;
			aDeltaM = monomers*delta_m;
			aDeltaP = polymers*delta_p;
			aBeta = polymers*monomers*beta;
			aB = polymerBonds*b;
			a = aLambda + aDeltaM + aBeta + aDeltaP + aB;
			waitTime = -1.0/(a*Math.log(r2));

			currentTime += waitTime;

			if(currentTime < endTime) {
				r1*=a;

				if(r1 < aLambda) {
					monomers += 1;
				} else if(r1 < aLambda + aDeltaM) {
					monomers-=1;
				} else if(r1 < aLambda + aDeltaM + aDeltaP) {				
					destroyedPolymer = samplePolymerLength(polymerLengths, polymers);
					polymerLengths.put(destroyedPolymer, Math.max(polymerLengths.get(destroyedPolymer)-1, 0));
					polymers-=1;
					polymerSubunits-=destroyedPolymer;
					polymerBonds-=destroyedPolymer-1;
				} else if(r1 < aLambda + aDeltaM + aDeltaP + aBeta) {		
					addedToPolymer = samplePolymerLength(polymerLengths, polymers);
					polymerLengths.put(addedToPolymer, Math.max(polymerLengths.get(addedToPolymer)-1, 0));
					
					if(polymerLengths.get(addedToPolymer+1) == null) {
						polymerLengths.put(addedToPolymer+1, 1);
					} else {
						polymerLengths.put(addedToPolymer+1, polymerLengths.get(addedToPolymer+1)+1);
					}

					monomers-=1;
					polymerBonds+=1;
					polymerSubunits+=1;
				} else {
					breakingPolymer = samplePolymerBond(polymerLengths, polymerBonds);
					
					polymerLengths.put(breakingPolymer, Math.max(polymerLengths.get(breakingPolymer)-1, 0));
					polymers-=1;
					polymerSubunits-=breakingPolymer;
					polymerBonds-=breakingPolymer-1;
					
					breakPosition = (int)((breakingPolymer-1)*Math.random())+1;
					largeFragment = Math.max(breakPosition, breakingPolymer-breakPosition);
					smallFragment = Math.min(breakPosition, breakingPolymer-breakPosition);

					if(largeFragment >= polymerThreshold) {
						polymers+=1;
						polymerBonds+=largeFragment-1;
						polymerSubunits+=largeFragment;
						
						if(polymerLengths.get(largeFragment) == null) {
							polymerLengths.put(largeFragment, 1);
						} else {
							polymerLengths.put(largeFragment, polymerLengths.get(largeFragment)+1);
						}
					} else {
						monomers += largeFragment;
					}

					if(smallFragment >= polymerThreshold) {
						polymers+=1;
						polymerBonds+=smallFragment-1;
						polymerSubunits+=smallFragment;
						
						if(polymerLengths.get(smallFragment) == null) {
							polymerLengths.put(smallFragment, 1);
						} else {
							polymerLengths.put(smallFragment, polymerLengths.get(smallFragment)+1);
						}
					} else {
						monomers += smallFragment;
					}
				}
			}
		}

		return polymerLengths;
	}

	public static int[][] prionGillespieInfectiousNoninfectious(double[] sampleTimes, double lambda, double delta_m, double beta, double delta_p, double b, int polymerThreshold, HashMap<Integer, Integer> polymerLengths0, int monomers0) {
		int[][] sampledValues = new int[sampleTimes.length][4];
		int counter = 0;
		
		// variables to keep track of state
		double currentTime = 0;
		HashMap<Integer, Integer> polymerLengths = (HashMap<Integer, Integer>)(polymerLengths0.clone());
		int polymers = 0;
		for(int i: polymerLengths.values())
			polymers+=i;
		int monomers = monomers0;
		int polymerSubunits = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
			polymerSubunits+=entry.key()*entry.value();
		int polymerBonds = 0;
		for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
			polymerBonds+=(entry.key()-1)*entry.value();


		// variables to keep track of probabilities
		double aLambda, aDeltaM, aBeta, aDeltaP, aB, a;
		double r1, r2, waitTime;
		int destroyedPolymer, addedToPolymer, breakingPolymer, breakPosition, largeFragment, smallFragment;

		while(counter < sampleTimes.length) {
			r1 = Math.random();
			r2 = Math.random();

			aLambda = lambda;
			aDeltaM = monomers*delta_m;
			aDeltaP = polymers*delta_p;
			aBeta = polymers*monomers*beta;
			aB = polymerBonds*b;
			a = aLambda + aDeltaM + aBeta + aDeltaP + aB;
			waitTime = -1.0/(a*Math.log(r2));
			
			while(currentTime < sampleTimes[counter] && currentTime + waitTime >= sampleTimes[counter]) {
				sampledValues[counter][0] = monomers;
				sampledValues[counter][1] = polymers;
				sampledValues[counter][2] = polymerSubunits;
				sampledValues[counter][3] = sampleTimes[counter];
				counter++;
  			}

  			currentTime += waitTime;

			if(currentTime < endTime) {
				r1*=a;

				if(r1 < aLambda) {
					monomers += 1;
				} else if(r1 < aLambda + aDeltaM) {
					monomers-=1;
				} else if(r1 < aLambda + aDeltaM + aDeltaP) {				
					destroyedPolymer = samplePolymerLength(polymerLengths, polymers);
					polymerLengths.put(destroyedPolymer, Math.max(polymerLengths.get(destroyedPolymer)-1, 0));
					polymers-=1;
					polymerSubunits-=destroyedPolymer;
					polymerBonds-=destroyedPolymer-1;
				} else if(r1 < aLambda + aDeltaM + aDeltaP + aBeta) {		
					addedToPolymer = samplePolymerLength(polymerLengths, polymers);
					polymerLengths.put(addedToPolymer, Math.max(polymerLengths.get(addedToPolymer)-1, 0));
					
					if(polymerLengths.get(addedToPolymer+1) == null) {
						polymerLengths.put(addedToPolymer+1, 1);
					} else {
						polymerLengths.put(addedToPolymer+1, polymerLengths.get(addedToPolymer+1)+1);
					}

					monomers-=1;
					polymerBonds+=1;
					polymerSubunits+=1;
				} else {
					breakingPolymer = samplePolymerBond(polymerLengths, polymerBonds);
					
					polymerLengths.put(breakingPolymer, Math.max(polymerLengths.get(breakingPolymer)-1, 0));
					polymers-=1;
					polymerSubunits-=breakingPolymer;
					polymerBonds-=breakingPolymer-1;
					
					breakPosition = (int)((breakingPolymer-1)*Math.random())+1;
					largeFragment = Math.max(breakPosition, breakingPolymer-breakPosition);
					smallFragment = Math.min(breakPosition, breakingPolymer-breakPosition);

					if(largeFragment >= polymerThreshold) {
						polymers+=1;
						polymerBonds+=largeFragment-1;
						polymerSubunits+=largeFragment;
						
						if(polymerLengths.get(largeFragment) == null) {
							polymerLengths.put(largeFragment, 1);
						} else {
							polymerLengths.put(largeFragment, polymerLengths.get(largeFragment)+1);
						}
					} else {
						monomers += largeFragment;
					}

					if(smallFragment >= polymerThreshold) {
						polymers+=1;
						polymerBonds+=smallFragment-1;
						polymerSubunits+=smallFragment;
						
						if(polymerLengths.get(smallFragment) == null) {
							polymerLengths.put(smallFragment, 1);
						} else {
							polymerLengths.put(smallFragment, polymerLengths.get(smallFragment)+1);
						}
					} else {
						monomers += smallFragment;
					}
				}
			}
		}

		return sampledValues;
	}

	public static int[][] prionGillespieInfectiousNonInfectious(double endTime, int numIterations, double lambda, double delta_m, double beta, double delta_p, double b, int polymerThreshold, HashMap<Integer, Integer> polymerLengths0, int monomers0) {
		
		int[][] sampledValues = new int[numIterations][3];

		// variables to keep track of state
		double currentTime;
		int polymers, monomers, polymerSubunits, polymerBonds;
		HashMap<Integer, Integer> polymerLengths = new HashMap<Integer, Integer>();

		// variables to keep track of probabilities
		double aLambda, aDeltaM, aBeta, aDeltaP, aB, a;
		double r1, r2, waitTime;
		int destroyedPolymer, addedToPolymer, breakingPolymer, breakPosition, largeFragment, smallFragment;

		for(int iteration=0; iteration<numIterations; iteration++) {

			currentTime = 0;

			polymerLengths.clear();
			for(Map.Entry<Integer, Integer> initialEntry: polymerLengths0.entrySet())
				polymerLengths.put(initialEntry.key(), initialEntry.value());

			polymers = 0;
			for(int i: polymerLengths.values())
				polymers+=i;

			polymerSubunits = 0;
			for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
				polymerSubunits+=entry.key()*entry.value();
			
			polymerBonds = 0;
			for(Map.Entry<Integer, Integer> entry: polymerLengths.entrySet())
				polymerBonds+=(entry.key()-1)*entry.value();

			monomers = monomers0;

			while(currentTime < endTime) {
				r1 = Math.random();
				r2 = Math.random();

				aLambda = lambda;
				aDeltaM = monomers*delta_m;
				aDeltaP = polymers*delta_p;
				aBeta = polymers*monomers*beta;
				aB = polymerBonds*b;
				a = aLambda + aDeltaM + aBeta + aDeltaP + aB;
				waitTime = -1.0/(a*Math.log(r2));

				currentTime += waitTime;

				if(currentTime < endTime) {
					r1*=a;

					if(r1 < aLambda) {
						monomers += 1;
					} else if(r1 < aLambda + aDeltaM) {
						monomers-=1;
					} else if(r1 < aLambda + aDeltaM + aDeltaP) {				
						destroyedPolymer = samplePolymerLength(polymerLengths, polymers);
						polymerLengths.put(destroyedPolymer, Math.max(polymerLengths.get(destroyedPolymer)-1, 0));
						polymers-=1;
						polymerSubunits-=destroyedPolymer;
						polymerBonds-=destroyedPolymer-1;
					} else if(r1 < aLambda + aDeltaM + aDeltaP + aBeta) {		
						addedToPolymer = samplePolymerLength(polymerLengths, polymers);
						polymerLengths.put(addedToPolymer, Math.max(polymerLengths.get(addedToPolymer)-1, 0));
						
						if(polymerLengths.get(addedToPolymer+1) == null) {
							polymerLengths.put(addedToPolymer+1, 1);
						} else {
							polymerLengths.put(addedToPolymer+1, polymerLengths.get(addedToPolymer+1)+1);
						}

						monomers-=1;
						polymerBonds+=1;
						polymerSubunits+=1;
					} else {
						breakingPolymer = samplePolymerBond(polymerLengths, polymerBonds);
						
						polymerLengths.put(breakingPolymer, Math.max(polymerLengths.get(breakingPolymer)-1, 0));
						polymers-=1;
						polymerSubunits-=breakingPolymer;
						polymerBonds-=breakingPolymer-1;
						
						breakPosition = (int)((breakingPolymer-1)*Math.random())+1;
						largeFragment = Math.max(breakPosition, breakingPolymer-breakPosition);
						smallFragment = Math.min(breakPosition, breakingPolymer-breakPosition);

						if(largeFragment >= polymerThreshold) {
							polymers+=1;
							polymerBonds+=largeFragment-1;
							polymerSubunits+=largeFragment;
							
							if(polymerLengths.get(largeFragment) == null) {
								polymerLengths.put(largeFragment, 1);
							} else {
								polymerLengths.put(largeFragment, polymerLengths.get(largeFragment)+1);
							}
						} else {
							monomers += largeFragment;
						}

						if(smallFragment >= polymerThreshold) {
							polymers+=1;
							polymerBonds+=smallFragment-1;
							polymerSubunits+=smallFragment;
							
							if(polymerLengths.get(smallFragment) == null) {
								polymerLengths.put(smallFragment, 1);
							} else {
								polymerLengths.put(smallFragment, polymerLengths.get(smallFragment)+1);
							}
						} else {
							monomers += smallFragment;
						}
					}
				}
			}

			sampledValues[iteration][0] = monomers;
			sampledValues[iteration][1] = polymers;
			sampledValues[iteration][2] = polymerSubunits;
 		}

		return sampledValues;
	}
}