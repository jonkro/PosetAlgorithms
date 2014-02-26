import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;


public class Poset 
{
	HashMap<Integer,Element> elements;
	HashSet<Integer> minima;
	HashSet<Integer> maxima;


	static int[][] fibonacciTriangle; 


	Poset()
	{
		elements = new HashMap<Integer,Element>();
		minima = new HashSet<Integer>();
		maxima = new HashSet<Integer>();
	}
	
	
	Poset(HashMap<Integer,Element> el, HashSet<Integer> min, HashSet<Integer> max)
	{
		elements = el;
		minima = min;
		maxima = max;
	}
	
	
	int size()
	{
		return elements.size();
	}
	
	
	boolean isEmpty()
	{
		return elements.isEmpty();
	}
	
	
	HashSet<Integer> getMinima()
	{
		return minima;
	}
	
	
	Element addElem(int i)
	{
		Element e = new Element(i);
		elements.put(i, e);
		minima.add(i);
		maxima.add(i);
		return e;
	}
	
	
	Element getElem(int i)
	{
		return elements.get(i);
	}
	
	
	void deleteElement(Element e)
	{
		elements.remove(e.getID());
		for(Integer i : e.getUpShadow())
		{
			elements.get(i).deleteSmaller(e.getID());
		}
		for(Integer i : e.getDownShadow())
		{
			elements.get(i).deleteGreater(e.getID());
		}
		update();
	}
	
	
	void deleteElement(int j)
	{
		deleteElement(elements.get(j));
	}
	
	
	void addSmallerRelation(Element s, Element g)
	{
		//System.out.println(g.isMinumum());
		if(g.isMinumum())
		{
			minima.remove(g.getID());
		}
		if(s.isMaximum())
		{
			maxima.remove(s.getID());
		}
		g.newSmaller(s.getID());
		s.newGreater(g.getID());
	}
	
	
	void addSmallerRelation(int i, int j)
	{
		if(!elements.containsKey(i)) addElem(i);
		if(!elements.containsKey(j)) addElem(j);
		Element s = elements.get(i);
		Element g = elements.get(j);
		addSmallerRelation(s,g);
	}
	
	
	boolean isSmaller(int i, int j)
	{
		boolean isSmaller = false;
		Vector<Integer> queue = new Vector<Integer>();
		queue.add(i);
		Integer current = null;
		while(!queue.isEmpty() && !isSmaller)
		{
			current = queue.remove(0);
			if(current.intValue() == j) isSmaller = true;
			for(int k : getElem(current).getUpShadow())
			{
				queue.add(k);
			}
		}
		return isSmaller;
	}
	
	
	void update()
	{
		minima.clear();
		maxima.clear();
		for(Element e : elements.values())
		{	
			if(e.isMinumum()) minima.add(e.getID());
			if(e.isMaximum()) maxima.add(e.getID());
		}
	}
	
	
	public Vector<LinearExtension> allLinearExtensions()
	{
		Vector<LinearExtension> ext = new Vector<LinearExtension>();
		LinearExtension next = new LinearExtension();
		extRec(next,this,ext);
		//System.out.println();
		return ext;
	}
	
	
	public Vector<LinearExtension> linExtKJumps(int k)
	{
		Vector<LinearExtension> allExt = allLinearExtensions();
		Vector<LinearExtension> kExt = new Vector<LinearExtension>();
		for(LinearExtension le : allExt)
		{
			if(le.getJumps() == k)
			{
				kExt.addElement(le);
			}
		}
		return kExt;
	}
	
	public Vector<LinearExtension> greedyLinExtKJumps(int k)
	{
		Vector<LinearExtension> allExt = allGreedyLinearExtensions();
		Vector<LinearExtension> kExt = new Vector<LinearExtension>();
		for(LinearExtension le : allExt)
		{
			if(le.getJumps() == k)
			{
				kExt.addElement(le);
			}
		}
		return kExt;
	}
	
	
	public String linExtJumpNr(int k)
	{
		Vector<LinearExtension> allLinExt = allLinearExtensions();
		StringBuffer str = new StringBuffer();
		for(LinearExtension le : allLinExt)
		{
			if(le.getJumps() == k)
			{
				//str.append(le.showMiddleChain());
				str.append(le.toString());
				str.append("\n");
			}
	
	
		}
		return str.toString();
	}
	
	
	public void extRec(LinearExtension prefix,Poset p,Vector<LinearExtension> ext)
	{
		LinearExtension next;
		//HashSet<Integer> help;
		Poset iter;
		if(!p.isEmpty())
		{
			for(Integer i : p.getMinima())
			{
				next = new LinearExtension(prefix);
				if(!prefix.isEmpty()) next.add(i,!isSmaller(prefix.lastElem(),i));
				else next.add(i);
				iter = p.copy();
				iter.deleteElement(i);
				//iter.update();
				extRec(next,iter,ext);
			}
		}
		else ext.add(prefix);
		//return next;
	}
	
	
	public Vector<LinearExtension> allGreedyLinearExtensions()
	{
		Vector<LinearExtension> ext = new Vector<LinearExtension>();
		LinearExtension next = new LinearExtension();
		greedyExtRec(next,this,ext);
		//System.out.println();
		return ext;
	}
	
	
	public void greedyExtRec(LinearExtension prefix,Poset p,Vector<LinearExtension> ext)
	{
		LinearExtension next;
		//HashSet<Integer> help;
		Poset iter;
		
		
		if(!p.isEmpty())
		{
		HashSet<Integer> options;
		if(prefix.isEmpty()||intersect(getElem(prefix.lastElem()).getUpShadow(), p.getMinima()).isEmpty())
		options = p.getMinima();
		else
		options = intersect(getElem(prefix.lastElem()).getUpShadow(), p.getMinima());
		for(Integer i : options)
		{
		next = new LinearExtension(prefix);
		if(!prefix.isEmpty()) next.add(i,!isSmaller(prefix.lastElem(),i));
		else next.add(i);
		iter = p.copy();
		iter.deleteElement(i);
		//iter.update();
		greedyExtRec(next,iter,ext);
		}
		}
		else ext.add(prefix);
		//return next;
	}
	
	
	public static HashSet<Integer> intersect(HashSet<Integer> a, HashSet<Integer> b)
	{
		HashSet<Integer> result = new HashSet<Integer>();
		for(int i : a)
		{
			if(b.contains(i))
			{
				result.add(i);
			}
		}
		return result;
	}
	
	
	public int countLinExt()
	{
		return iterLinExtCount(this,0);
	}
	
	
	int iterLinExtCount(Poset p, int k)
	{
		Poset iter;
		if(!p.isEmpty())
		{
			int r = 0;
			for(Integer i : p.getMinima())
			{
				iter = p.copy();
				iter.deleteElement(i);
				//iter.update();
				r += iterLinExtCount(iter,k);
			}
			return r;
		}
		else return 1;
	}
	
	
	Poset copy()
	{
	
	
		HashMap<Integer,Element> el = new HashMap<Integer,Element>();
		for(Integer i : elements.keySet())
		{
			el.put(i, elements.get(i).copy());
		}
		HashSet<Integer> min = new HashSet<Integer>();
		for(int i : minima)
		{
			min.add(i);
		}
		HashSet<Integer> max = new HashSet<Integer>();
		for(int i : maxima)
		{
			max.add(i);
		}
		return new Poset(el,min,max);
	}
	
	
	Poset copy(int offset)
	{
		HashMap<Integer,Element> el = new HashMap<Integer,Element>();
		for(Integer i : elements.keySet())
		{
			el.put(i+offset, elements.get(i).copy());
		}
		HashSet<Integer> min = new HashSet<Integer>();
		for(int i : minima)
		{
			min.add(i+offset);
		}
		HashSet<Integer> max = new HashSet<Integer>();
		for(int i : maxima)
		{
			max.add(i+offset);
		}
		return new Poset(el,min,max);
	}
	
	
	Poset parallelize(Poset q)
	{
		Poset r = this.copy();
		int offset = this.size();
		for(Integer i : q.elements.keySet())
		{
			r.elements.put(i+offset, q.elements.get(i).copy(offset));
		}
		for(int i : q.minima)
		{
			r.minima.add(i+offset);
		}
		for(int i : q.maxima)
		{
			r.maxima.add(i+offset);
		}
		return r;
	}
	
	
	static Poset chain(int r)
	{
		Poset c = new Poset();
		for(int i=1;i<r;i++)
		{
			c.addSmallerRelation(i, i+1);
		}
		return c;
	}
	
	
	public static String printExt(Vector<LinearExtension> extensions)
	{
		StringBuffer str = new StringBuffer();
		for(LinearExtension ext : extensions)
		{
			for(int i=0;i<ext.size()-1;i++)
			{
				str.append(ext.getElemAt(i)+",");
			}
			str.append(ext.lastElem()+"\n");
		}
		return str.toString();
	}
	
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		for(Element e : elements.values())
		{
			for(Integer i : e.getUpShadow())
			{
				str.append("("+e.getID()+","+i+"),");
			}
		}
		return str.toString();
	}
	
	
	public static Poset gridPoset(int n, int k)
	{
		Poset p = new Poset();
		for(int i=1;i<=n;i++)
		{
			for(int j=0;j<k;j++)
			{
				if(i<n) p.addSmallerRelation(i+n*j, i+n*j+1);
				if(j<k-1)
					p.addSmallerRelation(i+n*j, i+n*(j+1));
			}
		}
		return p;
	}
	
	
	public static Poset ladder(int x, int y)
	{
		Poset p = new Poset();
		for(int i=1;i<y;i++)
		{
			p.addSmallerRelation(i, i+1);
			p.addSmallerRelation(x+i, x+i+1);
			p.addSmallerRelation(i, x+i);
		}
		p.addSmallerRelation(y, x+y);
		for(int i=y;i<x;i++)
		{
			p.addSmallerRelation(i, i+1);
		}
		return p;
	}
	
	
	public static Poset ladder(int n)
	{
		return ladder(n,n);
	}
	
	
	public static Poset twoChains(int x, int y)
	{
		Poset p = new Poset();
		for(int i=1;i<y;i++)
		{
			p.addSmallerRelation(i, i+1);
			p.addSmallerRelation(x+i, x+i+1);
			//p.addSmallerRelation(i, x+i);
		}
		//p.addSmallerRelation(y, x+y);
		for(int i=y;i<x;i++)
		{
			p.addSmallerRelation(i, i+1);
		}
		return p;
	}
	
	
	public static Poset fibonacciLadder(int n)
	{
		int k = (int) Math.ceil(((double)n)/2);
		Poset p = new Poset();
		if(n>0) p.addElem(1);
		if(n>1) p.addElem(1+k);
		for(int i=1;i<k;i++)
		{
			p.addSmallerRelation(i, i+1);
			if(i+k+1<=n)
			{
				p.addSmallerRelation(i+k, i+k+1);
				p.addSmallerRelation(i, i+k+1);
			}
			if(i>=2)
			{
				p.addSmallerRelation(i+k-1,i+1);
			}
		}
		return p;
	}
	
	
	public static Poset completeBinTree(int n)
	{
		Poset binTree = new Poset();
		binTree.addElem(1);
		for(int i=1; i<n; i++)
		{
			for(int j=(int) Math.pow(2,i); j<Math.pow(2,i+1); j++)
			binTree.addSmallerRelation(j/2, j);
		}
		return binTree;
	}
	
	
	public static Poset oneSidedCompleteBinTree(int n)
	{
		Poset binTree = new Poset();
		binTree.addElem(1);
		for(int i=1; i<n; i++)
		{
			for(int j=2*i; j<2*i+2; j++)
			binTree.addSmallerRelation(2*i-1, j);
		}
		return binTree;
	}
	
	
	static int gridLinExtCalc(int l, int r)
	{
		int[] numExp = new int[r*l];
		for(int i=0;i<r*l;i++)
		{
			numExp[i] = 1;
		}
		int[] denomExp = new int[l+r-1];
		for(int i=0;i<l+r-1;i++)
		{
			denomExp[i] = Math.min(Math.min(i+1,r+l-i-1), Math.min(r,l)); 
		}
		Fraction res = new Fraction(numExp,denomExp);
		return res.value();
	}
	
	
	public static int dDimCatalan(int d, int n)
	{
		int[] numExp = new int[d*n];
		for(int i=0;i<d*n;i++)
		{
			numExp[i] = Math.max(1, d-i);
		}
		int[] denomExp = new int[d+n-1];
		for(int i=0;i<d+n-1;i++)
		{
			denomExp[i] = Math.min(d, d+n-1-i);
		}
		Fraction res = new Fraction(numExp,denomExp);
		return res.value();
	}
	
	
	public static int[] jumpCount(Vector<LinearExtension> extensions)
	{
		int[] jumpNumbers = new int[extensions.get(0).size()];
		for(LinearExtension ext : extensions)
		{
			jumpNumbers[ext.getJumps()]++;
		}
		return jumpNumbers;
	}
	
	
	
	public static String toString(int[] list)
	{
		StringBuffer str = new StringBuffer();
		for(int i=0; i<list.length; i++)
		str.append(list[i] + "\t");
		return str.toString();
	}
	
	
	public static String ladderJumpTriangle(int s)
	{
		StringBuffer str = new StringBuffer();
		for(int n = 1; n <= s; n++)
		{
			for(int k = 0; k<= 2*(n-1); k++)
			{
				str.append(a(n,k));
				str.append("\t");
			}
			str.append("\n");
		}
		return str.toString();
	}
	
	
	public static String fibonacciJumpTriangle(int s)
	{
		StringBuffer str = new StringBuffer();
		for(int n = 1; n <= s; n++)
		{
			str.append(fibonacciJumps(n));
			str.append("\n");
		}
		return str.toString();
	}
	
	
	public static int[][] fibonacciJumpArray(int s)
	{
		int[][] triangle = new int[s][s];
		int[] line = null;
		for(int n = 1; n <= s; n++)
		{
			Poset p = fibonacciLadder(n);
			Vector<LinearExtension> linExt = p.allLinearExtensions();
			line = jumpCount(linExt);
			for(int k=0; k<n; k++)
			{
				triangle[n-1][k] = line[k];
			}
		}
		fibonacciTriangle = triangle;
		return triangle;
	}
	
	
	public static int fibonacciJumps(int n, int k)
	{
		if(n<=0 || k < 0 || k>=n)
			return 0;
		else if(fibonacciTriangle==null || k>fibonacciTriangle[0].length-1 || n>fibonacciTriangle[0].length)
		{
			fibonacciJumpArray(Math.max(n, k));
		}
		return fibonacciTriangle[n-1][k];
	}
	
	
	public static String ladderJumpSymmetricTriangle(int s)
	{
		StringBuffer str = new StringBuffer();
		for(int n = 1; n <= s; n++)
		{
			for(int i=0;i<s-n+1;i++)
				str.append("\t");
			for(int k = 0; k<= 2*(n-1); k++)
			{
				str.append(a(n,k));
				str.append("\t");
			}
			str.append("\n");
		}
		return str.toString();
	}
	
	
	private static int h(int n, int j, int i)
	{
		int res = 0;
		if(i==n&&j==0) res = 1;
		else if(i+j<n && j>0)
		{
			int m = Math.min((i+j),(n-j));
			for(int h=i+1;h<=m;h++)
			{
				res += Fraction.simpleFrac(h-i,n-h).multiply(Fraction.binomial(h,i)).multiply(Fraction.binomial(n-h, j)).multiply(Fraction.binomial(n-h, i+j-h)).value();
			}
		}
		return res;
	}
	
	
	public static int a(int n, int k)
	{
		int a = 0;
		int j = 0;
		for(int i=0;i<=n;i++)
		{
			if((k+1-i)%2 == 0) 
			{
				j = (k+1-i)/2;
				a += h(n,j,i);
			}
		}
		return a;
	}
	
	
	public static String diag(int n,int k, int l)
	{
		StringBuffer str = new StringBuffer();
		for(int i=0;i<l;i++)
		{
			str.append(a(n+i,k));
			str.append(",");
		}
		return str.toString();
	}
	
	
	public static String fibSeries(int k)
	{
		StringBuffer str = new StringBuffer();
		for(int i=0;i<=k;i++)
		{
		
			str.append(fibonacciLadder(i).allLinearExtensions().size());
			str.append(",");
		}
		return str.toString();
	}
	
	
	public static String fibonacciJumps(int k)
	{
		Poset p = fibonacciLadder(k);
		Vector<LinearExtension> linExt = p.allLinearExtensions();
		return toString(jumpCount(linExt));
	}
	
	
	public static int testRecurrence(int n, int k)
	{
		int result = fibonacciJumps(n-1,k-1);
		if(2*k-1==n) result--;
		int i = 2;
		while(k-i >= 0 && n-2*i >= 1)
		{
			result +=fibonacciJumps(n-2*i,k-i);
			if(2*k-1==n) result--;
			i++;
		}
		return result;
	}
	
	
	public static int recurrence2(int n, int k)
	{
		int result =0 ;
		if(n<=-1 || k<=-1)
			result = 0;
		else if((k==n-1 && k!=1) ||(k==0&&n==0))
			result = 1;
		else if (n==2 && k==1)
			result = 2;
		else if (n>=1 && k>=0 && k<=n-2) 
		{
			/*for(int b=0; b<=k;b++)
			{
				result += fibonacciJumps(n-2*b,k-b);
			}*/
			for (int a = 1; a <= k ; a++) 
			{
				result += recurrence2(n - a -2, k - a );
			}
			result += recurrence2(n - 2, k - 1);
		}
		return result;
	}
	
	
	static int twoChainJumps(int n, int m, int k)
	{
		int result = 0;
		if(n<=0||m<=0) result = 1;
		else
		{
			if(k==1) result = 2;
			else if(k==2) result = m+n-2;
			else if(k%2==1) 
				for(int l=0;l<=m-1-(k-1)/2;l++)
				result += 2*Fraction.binomial(((k-1)/2)-1+l,((k-1)/2)-1).value()*Fraction.binomial(n-1,(k-1)/2).value();
				else
				{
					result = Fraction.binomial(m-1,(k/2)-1).value()*Fraction.binomial(n-1,k/2).value();
					if(m-k/2>=1)
						for(int l=0;l<m-k/2;l++)
							result+=(m-(k/2)-l)*Fraction.binomial((k/2)+l-2, (k/2)-2).value()*Fraction.binomial(n-1, (k/2)-1).value();
				}
		}
		return result;
	}
	
	
	static int altTwoChainJumps(int n, int m, int k)
	{
		int result = 0;
		if(k%2==0)
		{
			result = Fraction.binValue(n-1, (k/2)-1)*Fraction.binValue(m-1, (k/2)) + Fraction.binValue(m-1, (k/2)-1)*Fraction.binValue(n-1, (k/2));
		}
		else
		{
			result = 2*Fraction.binValue(n-1, (k-1)/2)*Fraction.binValue(m-1, (k-1)/2);
		}
		return result;
	}
	
	
	
	static int twoChainGLEedges(int n, int m)
	{
	int result = 0;
	for(int i=1;i<=n;i++)
	for(int j=1;j<=m;j++)
	result += Fraction.binomial(i+j-2, i-1).value()*Fraction.binomial(n+m-i-j,n-i).value();
	return result;
	}
	
	
	static int twoChainGLEedges2(int n, int m)
	{
	int result = m;
	for(int k=1;k<=m;k++)
	for(int r=0;r<=k-1;r++)
	result += ((m-k+1)*(r+2)-1)*Fraction.binomial(k-1, r).value()*Fraction.binomial(n-1, r+1).value();
	return result;
	}
	
	
	static int chainMergeTrial(int r, Poset p, int k)
	{
	int result = 0;
	int sum = 0;
	int[] jumpVector = jumpCount(p.allLinearExtensions());
	for(int j=1; j<=r; j++)
	{
	sum = 0;
	for(int i=0;i<=j;i++)
	{
	if(k-2*j+i>=0 && k-2*j+i<p.size())
	sum += Fraction.binomial(k+2-2*j+i, i).value()*Fraction.binomial(p.size()-k+2*j-i-1, j-i).value()*jumpVector[k-2*j+i];
	}
	result += Fraction.binomial(r-1, j-1).value() * sum;
	}
	return result;
	}
	
	/*
	static int parallelPosets(Poset p, Poset q, int k)
	{
	int result = 0;
	int innerSum = 0;
	int[] pVector = jumpCount(p.allLinearExtensions());
	int[] qVector = jumpCount(q.allLinearExtensions());
	for(int l=0;l<q.size();l++)
	{
	for(int j=0;j<q.size();j++)
	{
	for(int m=0;m<=j;m++)
	{
	innerSum = 0;
	for(int i=0;i<=j+1;i++)
	{
	if(k-(l-m)-2*(j+1)+i >=0 &&k-(l-m)-2*(j+1)+i<=p.size()-1)
	innerSum += Fraction.binValue(k-(l-m)-2*j+i, i) * Fraction.binValue( p.size()-k+(l-m)+2*j-i+1, j+1-i) * pVector[k-(l-m)-2*(j+1)+i];
	}
	result += Fraction.binValue(l, m)*Fraction.binValue(q.size()-1-l, j-m)*qVector[l]*innerSum;
	}
	}
	}
	return result;
	}
	*/
	
	
	
	static int parallelPosets(Poset p, Poset q, int k)
	{
	return parallelPosets(jumpCount(p.allLinearExtensions()),jumpCount(q.allLinearExtensions()),k,p.size(),q.size());
	}
	
	
	
	static int parallelPosets(int[] p, int[] q, int k, int np, int nq)
	{
	int result = 0;
	int innerSum = 0;
	for(int l=0;l<nq;l++)
	{
	for(int j=0;j<nq;j++)
	{
	for(int m=0;m<=j;m++)
	{
	innerSum = 0;
	for(int i=0;i<=j+1;i++)
	{
	if(k-(l-m)-2*(j+1)+i >=0 &&k-(l-m)-2*(j+1)+i<=np-1)
	innerSum += Fraction.binValue(k-(l-m)-2*j+i, i) * Fraction.binValue( np-k+(l-m)+2*j-i+1, j+1-i) * p[k-(l-m)-2*(j+1)+i];
	}
	result += Fraction.binValue(l, m)*Fraction.binValue(nq-1-l, j-m)*q[l]*innerSum;
	}
	}
	}
	return result;
	}
	
	
	static int parallelPosets2(Poset p, Poset q, int k)
	{
	return parallelPosets2(jumpCount(p.allLinearExtensions()),jumpCount(q.allLinearExtensions()),k,p.size(),q.size());
	}
	
	
	static int parallelPosets3(Poset p, Poset q, int k)
	{
	return parallelPosets3(jumpCount(p.allLinearExtensions()),jumpCount(q.allLinearExtensions()),k,p.size(),q.size());
	}
	
	
	static int parallelPosets2(int[] p, int[] q, int k, int np, int nq)
	{
	int result = 0;
	int help = 0;
	for(int pCount=0;pCount<np;pCount++)
	{
	for(int qCount=0;qCount<nq;qCount++)
	{
	for(int i=0;i<=pCount;i++)
	{
	for(int j=0;j<=qCount;j++)
	{
	if(k-pCount-qCount%2==0)
	help = Fraction.binValue(np-1-pCount,((k-pCount-qCount)/2)-1) * Fraction.binValue(nq-1-qCount,(k-pCount-qCount)/2) + Fraction.binValue(np-1-pCount,(k-pCount-qCount)/2) * Fraction.binValue(nq-1-qCount,((k-pCount-qCount)/2)-1);
	else
	help = 2 * Fraction.binValue(np-1-pCount,(k-pCount-qCount-1)/2) * Fraction.binValue(nq-1-qCount,(k-pCount-qCount-1)/2); 
	result += Fraction.binValue(pCount, i) * Fraction.binValue(qCount, j) * help * p[pCount]*q[qCount];
	}
	}
	}
	}
	return result;
	}
	
	
	static int parallelPosets3(int[] p, int[] q, int k, int np, int nq)
	{
		int result = 0;
		int help = 0;
		for(int pCount=0; pCount<np; pCount++)
		{
			for(int qCount=0; qCount<nq; qCount++)
			{
				for(int b=0; b<=pCount; b++)
				{
					for(int r=0; r<=np-1-pCount; r++)
					{
						//help = Fraction.binValue(np-1-pCount,r)*Fraction.binValue(pCount, b)*Fraction.binValue(nq-1-qCount, k - pCount - qCount - 1 - r);
						result += Fraction.binValue(np-1-pCount,r)*Fraction.binValue(pCount, b)*Fraction.binValue(nq-1-qCount, k - pCount - qCount - 1 - r) * Fraction.binValue(qCount + 2, 2*r+b-k+pCount+qCount+2)*p[pCount]*q[qCount];
					}
					//result += Fraction.binValue(np-1-pCount,r)*Fraction.binValue(pCount, b)*Fraction.binValue(nq-1-qCount, k - pCount - qCount - 1 - r) * Fraction.binValue(qCount, 2*r+b-k+pCount+qCount+1)*p[pCount]*q[qCount];
				}
			}
		}
		return result;
	}
	

	
	public static Poset fencePoset(int n)
	{
		Poset p = new Poset();
		for(int i=1;i<n;i++)
		{
			if(i%2==1)
				p.addSmallerRelation(i, i+1);
			else
				p.addSmallerRelation(i+1, i);
				//p.addSmallerRelation(i, x+i);
		}
		return p;
	}
	
	
	public static void fenceTriangle(int max)
	{
		Poset fence = new Poset();
		fence.addElem(1);
		System.out.println(toString(jumpCount(fence.allLinearExtensions())));
		for(int n=2;n<=max;n++)
		{
			if(n%2==0)
				fence.addSmallerRelation(n-1, n);
			else
				fence.addSmallerRelation(n, n-1);
			System.out.println(toString(jumpCount(fence.allLinearExtensions())));
		}
	}
	
	
	// only works for EVEN max
	public static String efficientFenceTriangle(int max)
	{
		StringBuffer str = new StringBuffer();
		int[] singleton = {1};
		int[] path = {1,0};
		int[] vPosetJumps = {0,2,0};
		HashMap<IntTuple,int[]> parallelFences = new HashMap<IntTuple,int[]>();
		Vector<int[]> smallerFences = new Vector<int[]>();
		smallerFences.add(singleton);
		smallerFences.add(path);
		smallerFences.add(vPosetJumps);
		if(parallelFences.size()<max)
		{
			fenceJumps(max,parallelFences,smallerFences);
		}
		int[] currentSizeJumps;
		for(int n=1; n<=max; n++)
		{
			currentSizeJumps = smallerFences.get(n-1);
			for(int k=0;k<n-1;k++)
			{
				str.append(currentSizeJumps[k] + ",");
			}
			str.append(currentSizeJumps[n-1]+"\n");
		}
		return str.toString();
	}
	
	public static String parallelFences(int sum)
	{
		if(sum%2==1) sum++;
		sum+=2;
		StringBuffer str = new StringBuffer();
		int[] singleton = {1};
		int[] path = {1,0};
		int[] vPosetJumps = {0,2,0};
		HashMap<IntTuple,int[]> parallelFences = new HashMap<IntTuple,int[]>();
		Vector<int[]> smallerFences = new Vector<int[]>();
		smallerFences.add(singleton);
		smallerFences.add(path);
		smallerFences.add(vPosetJumps);
		if(parallelFences.size()<sum)
		{
			fenceJumps(sum,parallelFences,smallerFences);
		}
		int[] currentSizeJumps;
		for(IntTuple tuple : parallelFences.keySet())
		{
			str.append(tuple + ": ");
			currentSizeJumps = parallelFences.get(tuple);
			for(int k=0;k<currentSizeJumps.length-1;k++)
			{
				str.append(currentSizeJumps[k] + ",");
			}
			str.append(currentSizeJumps[currentSizeJumps.length-1]+"\n");
		}
		return str.toString();
	}
	
	public static int[] fenceJumps(int n)
	{
		int[] singleton = {1};
		int[] path = {1,0};
		int[] vPosetJumps = {0,2,0};
		HashMap<IntTuple,int[]> parallelFences = new HashMap<IntTuple,int[]>();
		Vector<int[]> smallerFences = new Vector<int[]>();
		smallerFences.add(singleton);
		smallerFences.add(path);
		smallerFences.add(vPosetJumps);
		return fenceJumps(n,parallelFences,smallerFences);
	}
	
	
	public static int[] fenceJumps(int n,HashMap<IntTuple,int[]> parallelFences, Vector<int[]> smallerFences)
	{
		//int[] jumps = null;
		int[] currentJumps = null;
		int[] knownJumps;
		IntTuple currentTuple;
		HashMap<Integer,HashMap<Integer,IntTuple>> tupleSet = new HashMap<Integer,HashMap<Integer,IntTuple>>();
		HashMap<Integer,IntTuple> currentMap = new HashMap<Integer,IntTuple>();
		int[] localP;
		int[] localQ;
		for(int i=2;i<=n-2;i++)
		{
			currentMap = new HashMap<Integer,IntTuple>();
			tupleSet.put(i, currentMap);
			for(int j=1;j<=n-2;j++)
			{
				currentMap.put(j, new IntTuple(i,j));
			}
			
		}	
		for(int localN=5;localN<=n;localN+=2)
		{
			//localN = i+2;
			currentJumps = new int[localN];
			
			
			int[] parallelJumps = null;
			for(int j=3;j<= localN-4; j+=2)
			{
				if(j <= localN-j-1)
					currentTuple = tupleSet.get(j).get(localN-j-1);
				else
					currentTuple = tupleSet.get(localN-j-1).get(j);
				if(!parallelFences.containsKey(currentTuple))
				{
					parallelJumps = new int[localN-1];
					for(int k2=0;k2<localN-1;k2++)
						parallelJumps[k2] = parallelPosets3(smallerFences.get(j-1),smallerFences.get(localN-j-2),k2,j,localN-j-1);
					parallelFences.put(currentTuple, parallelJumps);
				}
			}
		
		
			knownJumps = smallerFences.get(localN-3);
			for(int k=1;k<localN;k++)
			{
		
				if(k-1 < localN -2)
				{
					currentJumps[k] = 2 *(knownJumps[k-1]);
				}
				if(k>=2 && k-2 < localN-2)
				{
					currentJumps[k] += 2 * (k - 1) * knownJumps[k-2];
				}
				if(k>=3 && k-3 < localN-2)
				{
					currentJumps[k] += 2 * (localN - k ) * knownJumps[k-3]; 
				}
				for(int j=3;j<= localN-4; j+=2)
				{
					if(j <= localN-j-1)
						currentTuple = tupleSet.get(j).get(localN-j-1);
					else
						currentTuple = tupleSet.get(localN-j-1).get(j);
					parallelJumps = parallelFences.get(currentTuple);
					/*if(localN!=parallelJumps.length)
					{
					System.out.println(localN + " - " + parallelJumps.length);
					}*/
					if(k<=parallelJumps.length)
						currentJumps[k] += parallelJumps[k-1];
				}
			}
			smallerFences.add(localN -2 , null);
			smallerFences.add(localN-1,currentJumps);
		}
		if(n%2==0)
		{
			for(int localN=4;localN<=n;localN+=2)
			{
				currentJumps = new int[localN];	
				int[] parallelJumps = null;
				for(int j=3;j<= localN-3; j+=2)
				{
					currentTuple = tupleSet.get(j).get(localN-j-1);
				
					if(!parallelFences.containsKey(currentTuple))
					{
						parallelJumps = new int[localN-1];
						localP = smallerFences.get(j-1);
						localQ = smallerFences.get(localN-j-2);
						for(int k2=0;k2<localN-1;k2++)
							parallelJumps[k2] = parallelPosets3(localP,localQ,k2,j,localN-j-1);
						parallelFences.put(currentTuple, parallelJumps);
					}	
				}	
				knownJumps = smallerFences.get(localN-3);
				for(int k=1;k<localN;k++)
				{
					if(k-1 < localN -2)
					{
						currentJumps[k] = (knownJumps[k-1]);
					}
					if(k>=2 && k-2 < localN-2)
					{
						currentJumps[k] += (k - 1) * knownJumps[k-2];
					}
					if(k>=3 && k-3 < localN-2)
					{
						currentJumps[k] += (localN - k ) * knownJumps[k-3]; 
					}
					currentJumps[k] += smallerFences.get(localN-2)[k-1];
					for(int j=3;j<= localN-3; j+=2)
					{
						currentTuple = tupleSet.get(j).get(localN-j-1);
			
						parallelJumps = parallelFences.get(currentTuple);
						/*if(localN!=parallelJumps.length)
						{
							System.out.println(localN + " - " + parallelJumps.length);
						}*/
						if(k<=parallelJumps.length)
							currentJumps[k] += parallelJumps[k-1];
					}
				}
				if(smallerFences.size()==localN-1)
					smallerFences.add(currentJumps);
				else
					smallerFences.setElementAt(currentJumps, localN-1);
			}

		}
		return smallerFences.get(n-1);	
	}
	
	
	public static int[][] completeBinTreeJumps(int n)
	{
		//int[] result = {1};
		int[][] binTreeTable = new int[n][(int)Math.pow(2,n)];
		binTreeTable[0][0]=1;
		if(n>=2) 
		{
			binTreeTable[1][1]=2;
			for(int i=2;i<n;i++)
			{
				for(int k=0;k<Math.pow(2,i+1);k++)
					binTreeTable[i][k] = parallelPosets(binTreeTable[i-1],binTreeTable[i-1],k,(int) Math.pow(2,i)-1,(int) Math.pow(2,i)-1);
			}
		}
		return binTreeTable;
	}
	
	
	
	public static int[][] oneSidedCompleteBinTreeJumps(int n)
	{
		//int[] result = {1};
		int[][] binTreeTable = new int[n][(int)Math.pow(2,n)];
		binTreeTable[0][0]=1;
		if(n>=2) 
		{
			binTreeTable[1][1]=2;
			for(int i=2;i<n;i++)
			{
				for(int k=0;k<Math.pow(2,i+1);k++)
					binTreeTable[i][k] = parallelPosets(binTreeTable[i-1],binTreeTable[0],k,(int) Math.pow(2,i)-1,(int) Math.pow(2,i)-1);
			}
		}
		return binTreeTable;
	}
	
	
	
	public static void main(String[] args)
	{
		int n = 20;
		//int[] fenceJumps = fenceJumps(n);
		//System.out.println(toString(fenceJumps)); 
		System.out.println(parallelFences(n));
		//Poset fence = fencePoset(n);
		//System.out.println(toString(jumpCount(fence.allLinearExtensions())));
		//Poset openFence = fencePoset(n);
		//openFence.deleteElement(4);
		//System.out.println(toString(jumpCount(openFence.allLinearExtensions())));
	}



}