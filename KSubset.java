import java.util.Vector;


public class KSubset 
{
	boolean[] subset;
	
	KSubset(int n, int k)
	{
		subset = new boolean[n];
		if(k>n) k=n;
		for(int i=0;i<n;i++)
		{
			if(i<k) subset[i] = true;
			else subset[i] = false;
		}
	}
	
	KSubset(int n)
	{
		subset = new boolean[n];
		for(int i=0;i<n;i++)
		{
			subset[i] = false;
		}
	}
	
	boolean[] getSet()
	{
		return subset;
	}
	
	void addElem(int i)
	{
		subset[i] = true;
	}
	
	boolean isElem(int i)
	{
		return subset[i];
	}
	
	
	int getN()
	{
		return subset.length;
	}
	
	int getK()
	{
		int k=0;
		for(int i=0;i<getN();i++)
		{
			if(isElem(i)) k++;
		}
		return k;
	}
	
	static KSubset revDoorUnrank(int r, int n, int k)
	{
		int x = n;
		KSubset currentSubset = new KSubset(n);
		for(int i=k;i>0;i--)
		{
			//x = n;
			while(Fraction.binValue(x, i)>r) x--;
			currentSubset.addElem(x);
			r = Fraction.binValue(x+1, i) - r - 1;
		}
		return currentSubset;
	}
	
	public static Vector<KSubset> allKSubsets(int n, int k)
	{
		Vector<KSubset> allSubsets = new Vector<KSubset>();
		KSubset currentSubset;
		if(k==0)
		{
			currentSubset = new KSubset(n,0);
			allSubsets.add(currentSubset);
		}
		else
			for(int r=0;r<Fraction.binValue(n, k);r++)
			{
				currentSubset = revDoorUnrank(r,n,k);
				allSubsets.add(currentSubset);
			}
		return allSubsets;
	}
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		int k = getK();
		int j=0;
		str.append("[");
		for(int i=0;i<getN();i++)
		{
			if(isElem(i))
			{
				str.append(i);
				j++;
				if(j<k)
					str.append(",");
					
			}
		}
		str.append("]");
		return str.toString();
	}
	
	public static void main(String[] args)
	{
		//KSubset s = unrank(1,7,3);
		for(KSubset current : allKSubsets(7,3))
			System.out.println(current.toString());
	}
}
