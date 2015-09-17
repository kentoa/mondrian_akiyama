
import java.util.ArrayList;

public class MondrianQI {
    final static int k = 3, QI = 13, SA = 12, column = QI + SA + 1,rowSize = 8333;
    static ArrayList<double[][]> partitionSet = new ArrayList<double[][]>(); //mondrianÇ≈ï™äÑÇµÇΩÇ†Ç∆ÇÃÉfÅ[É^Çäiî[Ç∑ÇÈ
    
    public static void main(String[] args){
	String inputFilepath="/Users/akiyama/Desktop/PWS/GIJI_2004zensho_s_dataset.csv";
	double[][] data = new double[rowSize][QI + SA]; //ì¸óÕÉfÅ[É^
	ArrayList<double[]> partition = new ArrayList<double[]>(); //ÉåÉRÅ[ÉhÇÃèWçá
	
	ReadFile.FileScan(data,inputFilepath);
	AddRowId.addRowID(data,partition); //çsî‘çÜÇâ¡Ç¶ÇÈ
	mondrian(partition);
	SetValue(partitionSet);
	CSVOutput.CSVwrite(partitionSet);
    }
    
    static void SetValue(ArrayList<double[][]> partitionSet){ //ë„ï\ílÇ…íuÇ´ä∑Ç¶ÇÈä÷êî(äeQIÇÃíÜÇ≈ç≈ïpèoÇÃÇ‡ÇÃÇégÇ¡ÇƒíuÇ´ä∑Ç¶ÇÈ) 7,8,9,15Ç…íçñ⁄Ç∑ÇÈ
	int partitionSize = partitionSet.size();
	double [][][] candidatemax = new double[partitionSize][QI][2]; //äeëÆê´Ç≈ÇÃèoåªïpìxÇÃç≈ëÂíl,ÉCÉìÉfÉbÉNÉX
	int freq = 0;

	for(int i = 0; i < partitionSize; i++){
	    for(int j = 1; j < QI+1 ; j++){ //QIÇÕç≈ëÂïpìxÇÃÇ‡ÇÃÇ≈íuÇ´ä∑Ç¶
		freq=0;
		for(int k = 0; k < partitionSet.get(i).length-1; k++){
		    if(partitionSet.get(i)[k][j] != partitionSet.get(i)[k+1][j]){
			if(candidatemax[i][j-1][0] < freq){
			    candidatemax[i][j-1][0] = freq;
			    candidatemax[i][j-1][1] = k;
			}
		    }
		    else freq++;
		}
	    }
	}
	
	for(int i = 0; i < partitionSize; i++){
	    for(int j = 1; j<QI + 1; j++){
		for(int k = 0; k < partitionSet.get(i).length; k++){
		    if(j < 7 || 9 < j)
			partitionSet.get(i)[k][j] = partitionSet.get(i)[(int)candidatemax[i][j-1][1]][j];
		}
	    }
	}
    }

    static void SwapSA(ArrayList<double[][]> partitionSet){
	double tmp;
	for(int i=0, partitionSetSize=partitionSet.size() ; i < partitionSetSize; i++){
	    for(int j =0; j<partitionSet.get(i).length-1; j++){
		if(partitionSet.get(i)[j][7] == partitionSet.get(i)[j+1][7] && partitionSet.get(i)[j][8] == partitionSet.get(i)[j+1][8] && partitionSet.get(i)[j][9] == partitionSet.get(i)[j+1][9]){
		    for(int k = QI+1 ; k < column; k++){
			tmp = partitionSet.get(i)[j+1][k];
			partitionSet.get(i)[j+1][k] = partitionSet.get(i)[j][k];
			partitionSet.get(i)[j+1][k] = tmp;
		    }
		}
	    }
	}
	
	for(int i=0, partitionSetSize=partitionSet.size() ; i < partitionSetSize; i++){
	    for(int j = QI+1; j< QI+SA+1; j++){ //SAÇÕïΩãœílÇ≈íuÇ´ä∑Ç¶
		double sum=0;
		if(j == 15 || j == 21){
		    for(int k = 0, num=partitionSet.get(i).length-1; k < num; k++){
			sum+=partitionSet.get(i)[k][j];
		    }
		    double average=sum/(double)(partitionSet.get(i).length-1);
		    for(int k = 0, num=partitionSet.get(i).length-1; k < num; k++){
			partitionSet.get(i)[k][j]=average;
		    }
		}
	    }
	}
    }
    static int mondrian(ArrayList<double[]> partition){ //çƒãAÇ≈ï™äÑÇ∑ÇÈ
	double splitVal;
	
	System.out.println("currentSIZE:" + partition.size());
	if(partition.size() < 2*k){
	    double[][] splitTable=new double[partition.size()][column];
	    for(int i = 0; i < partition.size(); i++){
		for(int j = 0;j < column; j++){
		    splitTable[i][j] = partition.get(i)[j];
		}
	    }
	    partitionSet.add(splitTable);
	    return 0;
	}else{
	    int dim = chooseDimension(partition); //ç≈ëÂílÇ∆ç≈è¨ílÇÃç∑Ç™ç≈ëÂÇ≈Ç†ÇÈëÆê´ÇëIÇ‘
	    System.out.println("dim:"+dim);
	    splitVal = findMedian(partition,dim); //íÜâõílÇíTÇ∑
	    
	    ArrayList<double[]> lhs = new ArrayList<double[]>(); //ëÆê´ÇÃílÇ™splitValÇÊÇËÇ‡è¨Ç≥Ç¢É^ÉvÉãÇäiî[
	    ArrayList<double[]> rhs = new ArrayList<double[]>(); //ëÆê´ÇÃílÇ™splitValÇÊÇËÇ‡ëÂÇ´Ç¢É^ÉvÉãÇäiî[
	    
	    for(int i = 0; i < partition.size(); i++){
		if(partition.get(i)[dim] <= splitVal) lhs.add(partition.get(i));
		else rhs.add(partition.get(i));
	    }
	    
	    System.out.println("lhs:rhs="+lhs.size()+":"+rhs.size());

	    if(partition.size()!=lhs.size()){
		mondrian(lhs);
	    }
	    else{ //partition.size()=>2*kÇ≈2*kå¬ÇÃ1Ç™Ç†ÇÈÇ∆Ç´í ÇÈ
		double[][] splitTable=new double[lhs.size()][column];
		for(int i=0;i<lhs.size();i++){
		    for(int j=0;j<column;j++)
			splitTable[i][j]=lhs.get(i)[j];
		}
		partitionSet.add(splitTable);
		return 0;
	    }
	    if(partition.size()!=rhs.size()){
		mondrian(rhs);
	    }else{
		double[][] splitTable=new double[rhs.size()][column];
		for(int i=0;i<rhs.size();i++){
		    for(int j=0;j<column;j++){
			splitTable[i][j]=rhs.get(i)[j];
		    }
		}
		partitionSet.add(splitTable);
		return 0;
	    }
	}
	return 0;
    }
    
    static int chooseDimension(ArrayList<double[]> partition){ //QI=7,8,9Ç…Ç¬Ç¢ÇƒÇµÇ©ï™äÑÇµÇ»Ç¢
	int partitionSize = partition.size();
	double[] index = new double[3];
	
	for(int i = 7; i < 10; i++){
	    Qsort(partition, i, 0, partition.size()-1);
	    for(int j = 0; j < partitionSize-1; j++){
		if(partition.get(j)[i] != partition.get(j+1)[i]){
		    index[i-7]++;
		}
	    }
	}
	return min(index) + 7; //äeëÆê´èWçáÇ…Ç†ÇÈà·Ç§ílÇéùÇ¬Ç‡ÇÃÇÃå¬êîÇ™ç≈ëÂÇ∆Ç»ÇÈÉCÉìÉfÉbÉNÉXÇï‘Ç∑
    }
    
    static int min(double[] diff){ //ç∑Ç™ç≈ëÂÇ∆Ç»ÇÈÇ‡ÇÃÇÃÉCÉìÉfÉbÉNÉXÇï‘Ç∑
	int index = -1, diffSize = diff.length;
	double min = Double.MAX_VALUE;
	for(int i = 0; i < diffSize; i++)
	    if(diff[i] < min){
		min = diff[i];
		index = i+1;
	    }
	return index;
    }
    
    static int max(double[] diff){ //ç∑Ç™ç≈ëÂÇ∆Ç»ÇÈÇ‡ÇÃÇÃÉCÉìÉfÉbÉNÉXÇï‘Ç∑
	int index = -1, diffSize = diff.length;
	double max = Double.MIN_VALUE;
	for(int i = 0; i < diffSize; i++)
	    if(max < diff[i]){
		max = diff[i];
		index = i+1;
	    }
	return index;
    }
    
    static double findMedian(ArrayList<double[]> partition,int dim){ //íÜâõílÇíTÇ∑
	double splitVal;
	int partitionSize = partition.size();
	Qsort(partition, dim, 0, partitionSize-1);
	
	if(partitionSize%2 == 0)
	    splitVal = (partition.get(partitionSize/2 -1)[dim] + partition.get(partitionSize/2)[dim])/2;
	else
	    splitVal = partition.get(partitionSize/2 +1)[dim];
	return splitVal;
    }
    
    static double findAverage(ArrayList<double[]> partition,int dim){ //ïΩãœílÇíTÇ∑
	double average, sum = 0;
	int partitionSize = partition.size();
	
	for(int i = 0; i < partitionSize; i++){
	    sum += partition.get(i)[dim];
	}
	
	average = sum/partitionSize;
	return average;
    }
    
    static void Qsort(ArrayList<double[]> partition,int dim,int left,int right){
	int i = left, j = right, k; //i:É\Å[ÉgÇ∑ÇÈîzóÒÇÃç≈è¨çsî‘çÜ,j:É\Å[ÉgÇ∑ÇÈëÆê´ÇÃç≈ëÂçsî‘çÜ
	double pivot, tmp;

	pivot=partition.get((left+right)/2)[dim]; //äÓèÄílÇîzóÒÇÃíÜâõïtãﬂÇ…Ç∆ÇÈ
	    
	while(true){
	    while(partition.get(i)[dim] < pivot) //pivotÇÊÇËëÂÇ´Ç¢ílÇ™èoÇÈÇ‹Ç≈iÇëùâ¡Ç≥ÇπÇÈ
		i++;
	    while(pivot < partition.get(j)[dim]) //pivotÇÊÇËè¨Ç≥Ç¢ílÇ™èoÇÈÇ‹Ç≈jÇå∏è≠Ç≥ÇπÇÈ
		j--;
	    if(j <= i) //j<=iÇ»ÇÁñ≥å¿ÉãÅ[ÉvÇ©ÇÁî≤ÇØÇÈ
		break;
	    for(k = 0; k < column; k++){ //swap
		tmp = partition.get(i)[k];
		partition.get(i)[k] = partition.get(j)[k];
		partition.get(j)[k] = tmp;
	    }
	    i++; //éüÇÃÉfÅ[É^
	    j--;
	}

	if(left<i-1) //äÓèÄílÇÃç∂Ç…2à»è„óvëfÇ™Ç†ÇÍÇŒç∂ÇÃîzóÒÇQÉ\Å[ÉgÇ∑ÇÈ
	    Qsort(partition, dim, left, i-1);
	if(j+1<right) //äÓèÄílÇÃâEÇ…2à»è„óvëfÇ™Ç†ÇÍÇŒâEÇÃîzóÒÇQÉ\Å[ÉgÇ∑ÇÈ
	    Qsort(partition, dim, j+1, right);
    }
    
    
}
