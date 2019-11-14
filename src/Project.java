import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class Project extends JPanel {
	//file reading variables
	private static final String FILE = System.getProperty("user.dir") + "\\src\\data\\Random100.tsp";
	private static BufferedReader reader;
	
	//structures for storing information
	static ArrayList<Coordinate> locations = new ArrayList<Coordinate>();
	static ArrayList<Integer> perm = new ArrayList<Integer>();
	static ArrayList<ArrayList<Coordinate>> generation = new ArrayList<ArrayList<Coordinate>>();
	static ArrayList<Coordinate> firstParent = new ArrayList<Coordinate>();
	static ArrayList<Coordinate> secondParent = new ArrayList<Coordinate>();
	static ArrayList<Edge> edges = new ArrayList<Edge>();
	
	//variables for holding population and generation information
	static int populationSize = 100; //change for experiments
	static int populationCount = 0;
	static float generationCosts[] = new float[populationSize];
	
	static float highestCost;
	static float secondHighestCost;
	
	//function for reading in file information
	public static void getFileContent(File file) throws IOException {
		FileInputStream stream = new FileInputStream(file);
		reader = new BufferedReader(new InputStreamReader(stream));
		
		String line = null;
		while (!(line = reader.readLine()).equals("NODE_COORD_SECTION")) {
			//ignore metadata
		}
		
		while ((line = reader.readLine()) != null) {
			//split each location line on spaces to separate information
			String[] components = line.split("\\s+");
			
			//store each location using the Coordinate class
			Coordinate coord = new Coordinate();
			coord.setIdentifier(Integer.parseInt(components[0]));
			coord.setXCoord(Float.parseFloat(components[1]));
			coord.setYCoord(Float.parseFloat(components[2]));
			
			//add each location to the locations ArrayList
			locations.add(coord);
		}
		
		reader.close();
	}
	
	//function for calculating costs of all members of the population in the current generation
	static void calcCost(ArrayList<ArrayList<Coordinate>> gen) {
		float currentCost = 0;
		Coordinate firstCoord;
		Coordinate secondCoord;
		float x1;
		float y1;
		float x2;
		float y2;
		float radicand;
		
		//cycles through all members of the population
		for(int j = 0; j < gen.size(); j++) {
			//cycles through each point within each member of the population
			for (int i = 0; i < (gen.get(j).size() - 1); i++) {
				//get two adjacent locations in the permutation, starting at the first location
				firstCoord = gen.get(j).get(i);
				secondCoord = gen.get(j).get(i + 1);
				
				x1 = firstCoord.getXCoord();
				y1 = firstCoord.getYCoord();
				x2 = secondCoord.getXCoord();
				y2 = secondCoord.getYCoord();
				
				//calculate the distance between the two adjacent locations
				radicand = ((x2 - x1)*(x2 - x1)) + ((y2 - y1)*(y2 - y1));
				//add that distance to the total distance of the permutation's path so far
				currentCost += Math.sqrt(radicand);
			}	
			//store the member's cost
			generationCosts[j] = currentCost;
			//reset currentCost once the members total cost has been calculated
			currentCost = 0;
		
		}

	}
	
	//function using Heap's Algorithm to generate permutations
	//source: http://www.geeksforgeeks.org/heaps-algorithm-for-generating-permutations/
	static void heapPermutation(Integer[] perm2, int size, int n, int first) {
		
		//if size == 1 permuation is finished, will only run when populationCount < populationSize so not all permutations are generated
	    if (size == 1 && populationCount < populationSize) {
	    	Integer [] perm3 = new Integer[perm2.length + 2];
	    	ArrayList<Coordinate> generationMember = new ArrayList<Coordinate>();
		   
		   //add first location to beginning of permuation
	       perm3[0] = first;
	       //add permutation found with Heap's algorithm in between start/end location
	       for (int x = 0; x < perm2.length; x++) {
	    	   perm3[x + 1] = perm2[x]; 
	       }
	       //add first location to end of permutation so path is a cycle
	       perm3[perm3.length - 1] = first;
	       
	       //create an ArrayList of Coordinate objects from the permutation to become a member of the population
	       for(int j = 0; j < perm3.length; j++) {
	    	   Coordinate coord = new Coordinate();
	    	   coord.setIdentifier(perm3[j]);
	    	   coord.setXCoord(locations.get(perm3[j] - 1).getXCoord());
	    	   coord.setYCoord(locations.get(perm3[j] - 1).getYCoord());
	    	   generationMember.add(coord);
	       }
	       
	       //add the population member to the generation list and increment the populationCount
	       generation.add(generationMember);
	       populationCount++;
	       
	    } else if (populationCount >= populationSize) {
	    	return;
	    }
	 
	    for (int i = 0; i < size; i++) {
	    	//recursive call to heap algorithm
	        heapPermutation(perm2, size-1, n, first);
	 
	        // if size is odd, swap first and last element, else swap ith and last
	        if (size % 2 == 1) {
	    	   int temp = perm2[0];
	           perm2[0] = perm2[size-1];
	           perm2[size-1] = temp;
	        } else {
	    	   int temp = perm2[i];
	           perm2[i] = perm2[size-1];
	           perm2[size-1] = temp;
	        }
	    }
	}
	
	//function to find parents for the crossover functions to use
	public static void findParents(boolean isRandom, ArrayList<ArrayList<Coordinate>> gen) {
		if(isRandom) {
			//picks two random parents
			int firstParentIndex = ThreadLocalRandom.current().nextInt(0, gen.size());
			int secondParentIndex = ThreadLocalRandom.current().nextInt(0, gen.size());
			firstParent = gen.get(firstParentIndex);
			secondParent = gen.get(secondParentIndex);
		} else {
			//picks two best members of the population to be parents
			float bestCost = -1;
			int bestCostIndex = 0;
			float secondBestCost = -1;
			int secondBestCostIndex = 0;
			
			//find best member
			for(int x = 0; x < generationCosts.length; x++) {
				if((generationCosts[x] < bestCost) || (bestCost == -1)) {
					bestCost = generationCosts[x];
					bestCostIndex = x;
				}
			}
			firstParent = gen.get(bestCostIndex);
			
			//find second best member
			for(int y = 0; y < generationCosts.length; y++) {
				if(((generationCosts[y] < secondBestCost) || (secondBestCost == -1)) && (y != bestCostIndex)) {
					secondBestCost = generationCosts[y];
					secondBestCostIndex = y;
				}
			}
			secondParent = gen.get(secondBestCostIndex);
		}
	}
	
	//function to apply mutations to the population
	public static void Mutate(boolean isRandom, boolean isSwapWorst, int index) {
		if(isRandom) {
			//swap two random cities
			swapRandom(generation.get(index), index);
		} else if (isSwapWorst) {
			//swap two points that contribute to two worst edge costs
			swapWorst(generation.get(index), index);
		} else {
			//swap two halves of entire path
			swapHalves(generation.get(index), index);
		}
	}
	
	//function to swap two halves of the path passed to it
	public static void swapHalves(ArrayList<Coordinate> member, int index) {
		int halfwayIndex = (member.size()/2) - 1;
		ArrayList<Coordinate> child = new ArrayList<Coordinate>();
		float originalCost;
		float newCost;
		
		//add first location to mutation
		Coordinate coord = new Coordinate();
		coord.setIdentifier(1);
		coord.setXCoord(locations.get(0).getXCoord());
		coord.setYCoord(locations.get(0).getYCoord());
		child.add(coord);
		
		//add second half of original member to first half of mutation
		for(int x = halfwayIndex; x < (member.size() - 1); x++) {
			Coordinate coord2 = new Coordinate();
			coord2.setIdentifier(member.get(x).getIdentifier());
			coord2.setXCoord(member.get(x).getXCoord());
			coord2.setYCoord(member.get(x).getYCoord());
			child.add(coord2);
		}
		
		//add first half of original member to second half of mutation
		for(int x = 1; x < halfwayIndex; x++) {
			Coordinate coord2 = new Coordinate();
			coord2.setIdentifier(member.get(x).getIdentifier());
			coord2.setXCoord(member.get(x).getXCoord());
			coord2.setYCoord(member.get(x).getYCoord());
			child.add(coord2);
		}
		
		//add first location to the end of the mutation to complete the cycle
		child.add(coord);
		
		//calculate cost of original member and its mutation
		originalCost = calcCostOfMember(generation.get(index));
		newCost = calcCostOfMember(child);
		
		//replace original member with mutation only if the mutation is an improvement
		if(newCost < originalCost) {
			generation.remove(index);
			generation.add(child);
		}
	}
	
	//function to swap two random locations in a population member
	public static void swapRandom(ArrayList<Coordinate> member, int index) {
		int firstSwapIndex = 0;
		Coordinate firstSwapCoord;
		int firstSwapCoordId;
		float firstSwapCoordX;
		float firstSwapCoordY;
		
		int secondSwapIndex = 0;
		Coordinate secondSwapCoord;
		int secondSwapCoordId;
		float secondSwapCoordX;
		float secondSwapCoordY;
		
		float originalCost;
		float newCost;
		
		//pick two random locations to swap (but not the first/last location)
		firstSwapIndex = ThreadLocalRandom.current().nextInt(1, member.size());
		secondSwapIndex = ThreadLocalRandom.current().nextInt(1, member.size());
		
		//get data for first location
		firstSwapCoord = member.get(firstSwapIndex);
		firstSwapCoordId = firstSwapCoord.getIdentifier();
		firstSwapCoordX = firstSwapCoord.getXCoord();
		firstSwapCoordY = firstSwapCoord.getYCoord();
		
		//get data for second location
		secondSwapCoord = member.get(secondSwapIndex);
		secondSwapCoordId = secondSwapCoord.getIdentifier();
		secondSwapCoordX = secondSwapCoord.getXCoord();
		secondSwapCoordY = secondSwapCoord.getYCoord();
		
		//calculate cost of original member
		originalCost = calcCostOfMember(generation.get(index));
		
		//swap two locations
		generation.get(index).get(firstSwapIndex).setIdentifier(secondSwapCoordId);
		generation.get(index).get(firstSwapIndex).setXCoord(secondSwapCoordX);
		generation.get(index).get(firstSwapIndex).setYCoord(secondSwapCoordY);
		
		generation.get(index).get(secondSwapIndex).setIdentifier(firstSwapCoordId);
		generation.get(index).get(secondSwapIndex).setXCoord(firstSwapCoordX);
		generation.get(index).get(secondSwapIndex).setYCoord(firstSwapCoordY);
		
		//calculate cost of mutation
		newCost = calcCostOfMember(generation.get(index));
		
		//swap two locations back if the mutation is not better than the original member
		if(originalCost < newCost) {
			generation.get(index).get(secondSwapIndex).setIdentifier(secondSwapCoordId);
			generation.get(index).get(secondSwapIndex).setXCoord(secondSwapCoordX);
			generation.get(index).get(secondSwapIndex).setYCoord(secondSwapCoordY);
			
			generation.get(index).get(firstSwapIndex).setIdentifier(firstSwapCoordId);
			generation.get(index).get(firstSwapIndex).setXCoord(firstSwapCoordX);
			generation.get(index).get(firstSwapIndex).setYCoord(firstSwapCoordY);
		}
	}
	
	//function to swap two locations that contribute to the two highest costing edges
	public static void swapWorst(ArrayList<Coordinate> member, int index) {
		float edgeCost;
		float worstEdgeCost = -1;
		int firstSwapIndex = 0;
		float secondWorstEdgeCost = -1;
		int secondSwapIndex = 0;
		Coordinate firstCoord;
		Coordinate secondCoord;
		float x1;
		float y1;
		float x2;
		float y2;
		float radicand;
		
		float originalCost;
		float newCost;
		
		Coordinate firstSwapCoord;
		int firstSwapCoordId;
		float firstSwapCoordX;
		float firstSwapCoordY;
		Coordinate secondSwapCoord;
		int secondSwapCoordId;
		float secondSwapCoordX;
		float secondSwapCoordY;
		
		//find location contributing to worst edge
		for(int x = 0; x < (member.size() - 2); x++) {
			//get two adjacent locations in the permutation, starting at the first location
			firstCoord = member.get(x);
			secondCoord = member.get(x + 1);
			
			x1 = firstCoord.getXCoord();
			y1 = firstCoord.getYCoord();
			x2 = secondCoord.getXCoord();
			y2 = secondCoord.getYCoord();
			
			//calculate the distance between the two adjacent locations
			radicand = ((x2 - x1)*(x2 - x1)) + ((y2 - y1)*(y2 - y1));
			//add that distance to the total distance of the permutation's path so far
			edgeCost = (float) Math.sqrt(radicand);
			
			//find the worst cost and save the index one of the edges endpoints is at
			if(worstEdgeCost < edgeCost || worstEdgeCost == -1) {
				worstEdgeCost = edgeCost;
				firstSwapIndex = x + 1;
			}
		}
		
		//find location contributing to second worst edge
		for(int y = 0; y < (member.size() - 2); y++) {
			//get two adjacent locations in the permutation, starting at the first location
			firstCoord = member.get(y);
			secondCoord = member.get(y + 1);
			
			x1 = firstCoord.getXCoord();
			y1 = firstCoord.getYCoord();
			x2 = secondCoord.getXCoord();
			y2 = secondCoord.getYCoord();
			
			//calculate the distance between the two adjacent locations
			radicand = ((x2 - x1)*(x2 - x1)) + ((y2 - y1)*(y2 - y1));
			//add that distance to the total distance of the permutation's path so far
			edgeCost = (float) Math.sqrt(radicand);
			
			//find the second worst cost and save the index one of the edges endpoints is at
			if((secondWorstEdgeCost < edgeCost || secondWorstEdgeCost == -1) && ((y + 1) != firstSwapIndex)) {
				secondWorstEdgeCost = edgeCost;
				secondSwapIndex = y + 1;
			}
		}
		
		//get data for first location
		firstSwapCoord = member.get(firstSwapIndex);
		firstSwapCoordId = firstSwapCoord.getIdentifier();
		firstSwapCoordX = firstSwapCoord.getXCoord();
		firstSwapCoordY = firstSwapCoord.getYCoord();
		
		//get data for second location
		secondSwapCoord = member.get(secondSwapIndex);
		secondSwapCoordId = secondSwapCoord.getIdentifier();
		secondSwapCoordX = secondSwapCoord.getXCoord();
		secondSwapCoordY = secondSwapCoord.getYCoord();
		
		//calculate cost of original member
		originalCost = calcCostOfMember(generation.get(index));
		
		//swap two locations
		generation.get(index).get(firstSwapIndex).setIdentifier(secondSwapCoordId);
		generation.get(index).get(firstSwapIndex).setXCoord(secondSwapCoordX);
		generation.get(index).get(firstSwapIndex).setYCoord(secondSwapCoordY);
		
		generation.get(index).get(secondSwapIndex).setIdentifier(firstSwapCoordId);
		generation.get(index).get(secondSwapIndex).setXCoord(firstSwapCoordX);
		generation.get(index).get(secondSwapIndex).setYCoord(firstSwapCoordY);
		
		//calculate cost of mutation
		newCost = calcCostOfMember(generation.get(index));
		
		//swap two locations back if the mutation is not better than the original member
		if(originalCost < newCost) {
			generation.get(index).get(secondSwapIndex).setIdentifier(secondSwapCoordId);
			generation.get(index).get(secondSwapIndex).setXCoord(secondSwapCoordX);
			generation.get(index).get(secondSwapIndex).setYCoord(secondSwapCoordY);
			
			generation.get(index).get(firstSwapIndex).setIdentifier(firstSwapCoordId);
			generation.get(index).get(firstSwapIndex).setXCoord(firstSwapCoordX);
			generation.get(index).get(firstSwapIndex).setYCoord(firstSwapCoordY);
		}
	}
	
	//function to pick a member to be mutated
	public static int findMemberForMutation(boolean isRandom) {
		if (isRandom) {
			//generate a random member index
			return ThreadLocalRandom.current().nextInt(0, populationSize);
		} else {
			//get the index of the worst member of the population
			return findWorstMemberId();
		}
	}
	
	//function to find the best member of the population
	public static int findBestMember() {
		float bestCost = -1;
		int bestCostIndex = 0;
		
		for(int x = 0; x < generationCosts.length; x++) {
			if((generationCosts[x] < bestCost) || (bestCost == -1)) {
				bestCost = generationCosts[x];
				bestCostIndex = x;
			}
		}
		
		return bestCostIndex;
	}
	
	//function to create a list of edges from the member passed to it
	public static void buildEdges(ArrayList<Coordinate> member) {
		for(int x = 0; x < (member.size() - 1); x++) {
			Edge edge = new Edge();
			edge.setFirstIdentifier(member.get(x).getIdentifier());
			edge.setSecondIdentifier(member.get(x + 1).getIdentifier());
			edges.add(edge);
		}
	}
	
	//function to plot points and edges
	//help from https://www.java-forums.org/new-java/7995-how-plot-graph-java-given-samples.html
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        int h = getHeight();

        //print edges
        g2.setPaint(Color.blue);
        for(int i = 0; i < edges.size(); i++) {
        	int firstId = edges.get(i).getFirstIdentifier();
        	int secondId = edges.get(i).getSecondIdentifier();
        	
        	float x1 = locations.get(firstId - 1).getXCoord()*4;
        	float y1 = h - (locations.get(firstId - 1).getYCoord()*4);
        	float x2 = locations.get(secondId - 1).getXCoord()*4;
        	float y2 = h - (locations.get(secondId - 1).getYCoord()*4);
        	
        	Shape line = new Line2D.Float(x1, y1, x2, y2);
        	g2.draw(line);
        }
        
        //print points
        g2.setPaint(Color.green);
        for(int i = 0; i < locations.size(); i++) {
            float x = locations.get(i).getXCoord()*4;
            float y = h - (locations.get(i).getYCoord()*4);
            
            g2.fill(new Ellipse2D.Double(x-2, y-2, 4, 4));
        }
        
    }
    
    //function to find the id of the worst population member
    public static int findWorstMemberId() {
    	float worstCost = -1;
		int worstCostIndex = 0;
		
		for(int x = 0; x < generationCosts.length; x++) {
			if((generationCosts[x] > worstCost) || (worstCost == -1)) {
				worstCost = generationCosts[x];
				worstCostIndex = x;
			}
		}
		
		return worstCostIndex;
    }
    
    //function to calculate the cost of one member
    static float calcCostOfMember(ArrayList<Coordinate> member) {
		float cost = 0;
		Coordinate firstCoord;
		Coordinate secondCoord;
		float x1;
		float y1;
		float x2;
		float y2;
		float radicand;
		
		for (int i = 0; i < (member.size() - 1); i++) {
			//get two adjacent locations in the permutation, starting at the first location
			firstCoord = member.get(i);
			secondCoord = member.get(i + 1);
			
			x1 = firstCoord.getXCoord();
			y1 = firstCoord.getYCoord();
			x2 = secondCoord.getXCoord();
			y2 = secondCoord.getYCoord();
			
			//calculate the distance between the two adjacent locations
			radicand = ((x2 - x1)*(x2 - x1)) + ((y2 - y1)*(y2 - y1));
			//add that distance to the total distance of the permutation's path so far
			cost += Math.sqrt(radicand);
		}	
		
		return cost;

	}
    
    //crossover function using one crossover point
    public static void onePointCrossover(ArrayList<Coordinate> first, ArrayList<Coordinate> second) {
    	//get random crossover point
    	int crossoverPoint = ThreadLocalRandom.current().nextInt(0, (first.size() - 2));
    	
    	int y = crossoverPoint + 1;
    	ArrayList<Coordinate> child = new ArrayList<Coordinate>();
    	boolean added[] = new boolean[first.size() - 1];
    	int unaddedIndex = 0;
    	int worstId = 0;
    	
    	//initialize all locations as not added
    	for(int j = 0; j < (first.size() - 1); j++) {
    		added[j] = false;
    	}
    	
    	//add all locations from the first parent to the child up until the crossover point
    	for(int x = 0; x <= crossoverPoint; x++) {
    		Coordinate coord = new Coordinate();
    		coord.setIdentifier(first.get(x).getIdentifier());
    		coord.setXCoord(first.get(x).getXCoord());
    		coord.setYCoord(first.get(x).getYCoord());
    		
    		child.add(coord);
    		added[first.get(x).getIdentifier() - 1] = true;
    	}
		
    	//complete the child
    	while(y <= (second.size() - 2)) {
    		
    		//add all locations from the second parent to the child after the crossover point if the locations have not been added yet
    		if(added[(second.get(y).getIdentifier() - 1)] == false) {
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(second.get(y).getIdentifier());
        		coord.setXCoord(second.get(y).getXCoord());
        		coord.setYCoord(second.get(y).getYCoord());
        		
        		child.add(coord);
        		added[second.get(y).getIdentifier() - 1] = true;
    		} else {    
    			//if the location from the second parent has already been added, add the next unadded location
    			for(int k = 0; k < added.length; k++) {
    				if(added[k] == false) {
    					unaddedIndex = k;
    					break;
    				}
    			}
    			
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(unaddedIndex + 1);
        		coord.setXCoord(locations.get(unaddedIndex).getXCoord());
        		coord.setYCoord(locations.get(unaddedIndex).getYCoord());
        		
        		child.add(coord);
        		added[unaddedIndex] = true;
    		}
    		y++;
    	}
    	
    	//add the first location to the end of the child to complete the loop
    	Coordinate coord = new Coordinate();
		coord.setIdentifier(1);
		coord.setXCoord(locations.get(0).getXCoord());
		coord.setYCoord(locations.get(0).getYCoord());
		child.add(coord);
		
		//find the worst member
		worstId = findWorstMemberId();
		//replace the worst member with the new child if the new child is better
    	if(calcCostOfMember(child) < calcCostOfMember(generation.get(worstId))) {
	    	generation.remove(worstId);
	    	generation.add(child);
    	}
    }
   
  //crossover function using two crossover points
    public static void twoPointCrossover(ArrayList<Coordinate> first, ArrayList<Coordinate> second) {
    	//get first random crossover point
    	int crossoverPoint1 = ThreadLocalRandom.current().nextInt(0, (first.size() - 3));
    	int y = crossoverPoint1 + 1;
    	//get second random crossover point
    	int crossoverPoint2 = ThreadLocalRandom.current().nextInt(y, (first.size() - 2));
    	int z = crossoverPoint2 + 1;
    	
    	ArrayList<Coordinate> child = new ArrayList<Coordinate>();
    	boolean added[] = new boolean[first.size() - 1];
    	int unaddedIndex = 0;
    	int worstId = 0;
    	
    	//initialize all locations as not added
    	for(int j = 0; j < (first.size() - 1); j++) {
    		added[j] = false;
    	}
    	
    	//add all locations from the first parent to the child up until the first crossover point
    	for(int x = 0; x <= crossoverPoint1; x++) {
    		Coordinate coord = new Coordinate();
    		coord.setIdentifier(first.get(x).getIdentifier());
    		coord.setXCoord(first.get(x).getXCoord());
    		coord.setYCoord(first.get(x).getYCoord());
    		
    		child.add(coord);
    		added[first.get(x).getIdentifier() - 1] = true;
    	}
		
    	//continue the child
    	while(y <= crossoverPoint2) {
    		//add all locations from the second parent to the child up until the second crossover point
    		if(added[(second.get(y).getIdentifier() - 1)] == false) {
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(second.get(y).getIdentifier());
        		coord.setXCoord(second.get(y).getXCoord());
        		coord.setYCoord(second.get(y).getYCoord());
        		
        		child.add(coord);
        		added[second.get(y).getIdentifier() - 1] = true;
    		} else {    		
    			//if the location from the second parent has already been added, add the next unadded location
    			for(int k = 0; k < added.length; k++) {
    				if(added[k] == false) {
    					unaddedIndex = k;
    					break;
    				}
    			}
    			
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(unaddedIndex + 1);
        		coord.setXCoord(locations.get(unaddedIndex).getXCoord());
        		coord.setYCoord(locations.get(unaddedIndex).getYCoord());
        		
        		child.add(coord);
        		added[unaddedIndex] = true;
    		}
    		y++;
    	}
    	
    	//complete the child
    	while(z <= (second.size() - 2)) {
    		//add all locations from the first parent to the child from the second crossover point until the end if the locations have not been added yet
    		if(added[(first.get(z).getIdentifier() - 1)] == false) {
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(first.get(z).getIdentifier());
        		coord.setXCoord(first.get(z).getXCoord());
        		coord.setYCoord(first.get(z).getYCoord());
        		
        		child.add(coord);
        		added[first.get(z).getIdentifier() - 1] = true;
    		} else {    
    			//if the location from the second parent has already been added, add the next unadded location
    			for(int k = 0; k < added.length; k++) {
    				if(added[k] == false) {
    					unaddedIndex = k;
    					break;
    				}
    			}
    			Coordinate coord = new Coordinate();
        		coord.setIdentifier(unaddedIndex + 1);
        		coord.setXCoord(locations.get(unaddedIndex).getXCoord());
        		coord.setYCoord(locations.get(unaddedIndex).getYCoord());
        		
        		child.add(coord);
        		added[unaddedIndex] = true;
    		}
    		z++;
    	}
    	
    	//add the first location to the end of the child to complete the loop
    	Coordinate coord = new Coordinate();
		coord.setIdentifier(1);
		coord.setXCoord(locations.get(0).getXCoord());
		coord.setYCoord(locations.get(0).getYCoord());
		child.add(coord);
		
		//find the worst member
		worstId = findWorstMemberId();
		//replace the worst member with the new child if the new child is better
    	if(calcCostOfMember(child) < calcCostOfMember(generation.get(worstId))) {
	    	generation.remove(worstId);
	    	generation.add(child);
    	}
    }
	
	public static void main(String[] args) throws IOException {	
		//variable to keep track of first location in cycle
		int firstLocation = 0;
		
		//variables to change conditions of GA
		boolean isMutationSelectionRandom; //use this?
		boolean isMutationOperatorRandom;
		boolean isParentSelectionRandom;
		boolean isOnePointCrossover;
		boolean isMutationSwapWorst;
		
		//variables to control attributes of GA
		int terminationCondition = 100;
		float mutationRate = (float) .25;
		
		//variables for keeping track of indices
		int mutationMemberIndex;
		int bestIndex;
		
		//variables for statistics
		float min;
		float max;
		float avg = 0;
		float stddev = 0;
		
		//variables for calculating the time each algorithm takes
		long startTime;
		long endTime;
		
		//read in file information
		File file = new File(FILE); //line to change to change file
		getFileContent(file);
		
		//create list of location identifiers to use for generating permutations
		for (int n = 0; n < locations.size(); n++) {
			if(n == 0) {
				//store first location but don't use when generating permutations
				firstLocation = locations.get(n).getIdentifier();
			} else {
				//store all other locations to be used for generating permutations
				perm.add(locations.get(n).getIdentifier());
			}
		}
		
		//generate initial population and calculate their costs
		Integer [] perm2 = perm.toArray(new Integer[perm.size()]);
		heapPermutation(perm2, perm2.length, perm2.length, firstLocation);		
		calcCost(generation);
		
		//variables to alter conditions of the GA
		isParentSelectionRandom = false;
		isOnePointCrossover = true; //change for experiments
		isMutationSelectionRandom = false;
		isMutationOperatorRandom = true;
		isMutationSwapWorst = false;
		
		startTime = System.currentTimeMillis();
		
		//reproduce and mutate until termination condition is met, each loop represents a generation
		for(int x = 0; x < terminationCondition; x++) {	
			
			//create children
			for(int g = 0; g < (populationSize/2); g++) {
				findParents(isParentSelectionRandom, generation);
				if(isOnePointCrossover) {
					onePointCrossover(firstParent, secondParent);
				} else {
					twoPointCrossover(firstParent, secondParent);
				}
				calcCost(generation);
			}			
			
			//mutate
			for(int d = 0; d < (populationSize*mutationRate); d++) {
				mutationMemberIndex = findMemberForMutation(isMutationSelectionRandom);
				Mutate(isMutationOperatorRandom, isMutationSwapWorst, mutationMemberIndex);
				calcCost(generation);
			}
			
			//collect data for statistics
			min = calcCostOfMember(generation.get(findBestMember()));
			max = calcCostOfMember(generation.get(findWorstMemberId()));
			avg = 0;
			stddev = 0;
			for(int u = 0; u < generationCosts.length; u++) {
				avg += generationCosts[u];
			}
			avg = avg/generationCosts.length;
			
			for(int v = 0; v < generationCosts.length; v++) {
				stddev += (generationCosts[v] - avg)*(generationCosts[v] - avg);
			}
			stddev = stddev/generationCosts.length;
			stddev = (float) Math.sqrt(stddev);
			
			//print out statistics
			//System.out.println("Generation #" + (x + 1) + " Min: " + min + " Max: " + max + " Avg: " + avg + " Standard Deviation: " + stddev);
			System.out.println((x + 1) + "," + min + "," + max + "," + avg + "," + stddev);

		}
		
		//after GA is done get best path
		bestIndex = findBestMember();
		//print order of best member
		for(int a = 0; a < generation.get(bestIndex).size(); a++) {
			System.out.print(generation.get(bestIndex).get(a).getIdentifier() + " ");
		}
		System.out.println();
		
		//build and display best solution
		buildEdges(generation.get(bestIndex));
		
		//graph set of points and edges
		JFrame f = new JFrame("Visualization");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.add(new Project());
		f.setSize(450, 450);
		f.getContentPane().setSize(400,400);
		f.setLocation(200,200);
		f.setVisible(true);
		
		//calculate and print out runtime
		endTime = System.currentTimeMillis();
		System.out.println("Run time: " + (endTime - startTime));
	}

}