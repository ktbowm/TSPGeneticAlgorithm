//class for storing edge information
public class Edge {
	//class variables
	int firstPointId;
	int secondPointId;
	
	//setters
	public void setFirstIdentifier(int id) {
		firstPointId = id;
	}
	
	public void setSecondIdentifier(int id) {
		secondPointId = id;
	}
	
	//getters
	public int getFirstIdentifier() {
		return firstPointId;
	}
	
	public int getSecondIdentifier() {
		return secondPointId;
	}
	
	//printing functions
	public void printFirstIdentifier() {
		System.out.println(firstPointId);
	}
	
	public void printSecondIdentifier() {
		System.out.println(secondPointId);
	}

}