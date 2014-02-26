import java.util.HashMap;
import java.util.Vector;


public class Permutation 
{
	int n;
	int[] permutation;
	int[] direction;
	
	public Permutation(int size)
	{
		n = size;
		permutation = new int[n];
		direction = new int[n];
		for(int i=0;i<n;i++)
		{
			permutation[i] = i;
			direction[i] = -1;
		}
	}
	
	public Permutation(int[] p)
	{
		n = p.length;
		permutation = p;
	}
	
	int permParity(int m)
	{
		int parity=0;
		for(int j=0;j<m-1;j++)
		{
			for(int l=j+1;j<m;l++)
			{
				if(permutation[j]>permutation[l])
					parity++;
			}
		}
		return parity;
	}
	
	public int getEntry(int i)
	{
		return permutation[i];
	}
	public void setEntry(int pos, int val)
	{
		permutation[pos] = val;
	}
	
	public int getN()
	{
		return n;
	}
	
	void swap(int pos1, int pos2)
	{
		int help = permutation[pos1];
		permutation[pos1] = permutation[pos2];
		permutation[pos2] = help;
	}
	
	void permTrotterJohnsonSuccessor()
	{
		int s = 0;
		boolean done = false;
		int m = n-1;
		int[] r = new int[n];
		int d = 0;
		int par = 0;
		for(int i=0;i<n;i++)
		{
			r[i] = permutation[i];
		}
		Permutation rho = new Permutation(r);
		while(m>0 && !done)
		{
			d = 0;
			while(rho.getEntry(d)!=m-1)
				d++;
			for(int i=d;d<m;d++)
				rho.setEntry(i, rho.getEntry(i+1));
			par = rho.permParity(m-1);
			if(par==1)
			{
				if(d==m) m--;
				else
				{
					swap(s+d,s+d+1);
					done = true;
				}
			}
			else 
			{
				if(d==0)
				{
					m--;
					s++;
				}
				else
				{
					swap(s+d,s+d+1);
					done = true;
				}
			}
		}
		if(m==0)
		{
			for(int i=0;i<n;i++)
			{
				permutation[i] = i;
			}
		}
	}
	
	void move()
	{
		//HashMap<Integer,Integer> positions = new HashMap<Integer,Integer>();
		int neighborIndex;
		int indexOfMaxMovableElement = -1;
		int maxMovableElement = -1;
		for(int i=0; i<n; i++)
		{
			neighborIndex = i+direction[permutation[i]];
			if(neighborIndex>=0 && neighborIndex<n)
			{
				if(permutation[i]>permutation[neighborIndex] && permutation[i]>maxMovableElement)
				{
					maxMovableElement = permutation[i];
					indexOfMaxMovableElement = i;
				}
			}
		}
		swap(indexOfMaxMovableElement,indexOfMaxMovableElement+direction[maxMovableElement]);
		for(int j=maxMovableElement+1; j<n; j++)
		{
			direction[j] *= -1;
		}
	}
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		str.append("[");
		for(int i=0;i<n-1;i++)
		{
			str.append(permutation[i]+",");
		}
		str.append(permutation[n-1]+"]");
		return str.toString();
	}
	
	static public String arrayToString(int[] array)
	{
		StringBuffer str = new StringBuffer();
		str.append("[");
		for(int i=0;i<array.length-1;i++)
		{
			str.append(array[i]+",");
		}
		str.append(array[array.length-1]+"]");
		return str.toString();
	}
	
	static public String arrayToString(boolean[] array)
	{
		StringBuffer str = new StringBuffer();
		str.append("[");
		for(int i=0;i<array.length-1;i++)
		{
			str.append(array[i]+",");
		}
		str.append(array[array.length-1]+"]");
		return str.toString();
	}
	
	
	static Vector<int[]> allPermutations(int n)
	{
		int nFactorial = 1;
		for(int i=2; i<=n; i++)			
			nFactorial *= i;
		int[] currentPerm;
		Vector<int[]> allPerm = new Vector<int[]>();
		Permutation p = new Permutation(n);

		for(int k=0;k<nFactorial;k++)
		{
			currentPerm = new int[n];
			if(k>0) p.move();
			for(int i=0; i<n; i++)
			{
				currentPerm[i] = p.getEntry(i);
			}
			allPerm.add(currentPerm);
		}
		return allPerm;
	}
	
	static Vector<int[]> allPermutationsWithKGaps(int n, int k)
	{
		return permutationsWithKGaps(allPermutations(n), k);
	}
	
	static Vector<int[]> permutationsWithKGaps(Vector<int[]> allPerms, int k)
	{
		Vector<int[]> permutationsWithKGaps = new Vector<int[]>();
		for(int[] perm : allPerms)
		{
			if(gapNr(perm) == k)
			{
				permutationsWithKGaps.add(perm);
			}
		}
		return permutationsWithKGaps;
	}
	
	static HashMap<boolean[],Vector<int[]>> permutationsBySignature(Vector<int[]> allPerms)
	{
		int n = allPerms.get(0).length;
		HashMap<boolean[],Vector<int[]>> permutationsBySignature = new HashMap<boolean[],Vector<int[]>>();
		Vector<KSubset> signSet;
		boolean[] currentSign;
		Vector<int[]> currentPerms;
		for(int k=0;k<=n-2;k++)
		{
			signSet = KSubset.allKSubsets(n-2, k);
			for(KSubset sign : signSet)
			{
				for(int[] perm : allPerms)
				{
					currentSign = sign.getSet();
					if(permHasSign(perm,currentSign))
					{
						//allPerms.remove(perm);
						if(!permutationsBySignature.containsKey(currentSign))
						{
							currentPerms = new Vector<int[]>();
							currentPerms.add(perm);
							permutationsBySignature.put(currentSign, currentPerms);
						}
						else
						{
							permutationsBySignature.get(currentSign).add(perm);
						}
					}
				}
			}
				
		}
		return permutationsBySignature;
	}
	
	public String getPermutationsBySign(HashMap<boolean[],Vector<int[]>> permutationsBySignature)
	{
		StringBuffer str = new StringBuffer();
		for(boolean[] sign : permutationsBySignature.keySet())
		{
			str.append(arrayToString(sign)+"\n");
			for(int[] perm : permutationsBySignature.get(sign))
			{
				str.append("\t"+arrayToString(perm)+"\n");
			}
		}
		return str.toString();
	}
	
	public static String noPermutationsBySign(HashMap<boolean[],Vector<int[]>> permutationsBySignature)
	{
		StringBuffer str = new StringBuffer();
		for(boolean[] sign : permutationsBySignature.keySet())
		{
			str.append(arrayToString(sign)+": "+permutationsBySignature.get(sign).size());
		}
		return str.toString();
	}
	
	public static void printNoPermutationsBySign(int n)
	{
		System.out.println(noPermutationsBySign(permutationsBySignature(allPermutations(n))));
	}
	
	
	static boolean permHasSign(int[] perm, boolean[] sign)
	{
		int[] positions = positions(perm);
		boolean permHasSign = true;
		for(int i=2;i<perm.length;i++)
		{
			//System.out.println(gapAtPos(i,perm,positions)!=sign[i-1]);
			if(gapAtPos(i,perm,positions)!=sign[i-2])
				permHasSign = false;
		}
		return permHasSign;
	}
	
	static int[] positions(int[] perm)
	{
		int[] pos = new int[perm.length];
		for(int i=0; i<perm.length; i++)
		{
			pos[perm[i]] = i;
		}
		return pos;
	}
	
	static boolean gapAtPos(int i, int[] perm, int[] positions)
	{
		boolean gapAtPos = false;
		if(perm[i]>0 && perm[i]<perm.length-1)
		{	
			gapAtPos = positions[perm[i]-1]<i && positions[perm[i]+1]<i;
		}
		return gapAtPos;
	}
	
	static int gapNr(int[] perm)
	{
		int gapNr = 0;
		int indexI = -1;
		int indexIs = -1;
		int indexIg = -1;
		for(int i=1;i<perm.length-1;i++)
		{
			for(int j=0; j<perm.length; j++)
			{
				if(perm[j] == i) indexI=j;
				else if(perm[j] == i-1) indexIs=j;
				else if(perm[j] == i+1) indexIg=j;
			}
			if(indexI>Math.max(indexIs,indexIg)) gapNr++;
		}
		return gapNr;
	}
	
	public static void main(String[] args)
	{
		printNoPermutationsBySign(4);
	}
}
