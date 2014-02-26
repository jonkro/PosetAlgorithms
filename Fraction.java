
public class Fraction 
{
	int[] numExp;
	int[] denomExp;
	
	Fraction(int[] n, int[] d)
	{
		numExp = n;
		denomExp = d;
		reduce();
	}
	
	
	
	Fraction(int n, int d)
	{
		numExp = new int[n];
		if(n>=1)
			numExp[n-1]++;
		denomExp = new int[d];
		denomExp[d-1]++;
		reduce();
	}
	
	int value()
	{
		int e = 0;
		if(numExp.length>=1 && denomExp.length>=1)
		{
			reduce();
			e = 1;
			int num = numExp.length;
			int den = denomExp.length;
			for(int i=0;i<num;i++)
			{
				e *= Math.pow(i+1, numExp[i]);
			}
			for(int i=0;i<den;i++)
			{
				e = (int) e /(int) Math.pow(i+1, denomExp[i]);
			}
		}
		return e;
	}
	
	void reduce()
	{
		if(numExp.length>=1 && denomExp.length>=1)
		{
			int min = 0;
			int num = numExp.length;
			int den = denomExp.length;
			for(int i=0;i<Math.min(num,den);i++)
			{
				min = Math.min(numExp[i],denomExp[i]);
				numExp[i] = numExp[i] - min;
				denomExp[i] = denomExp[i] - min;
			}
			int j = num-1;
			int k;
			while(j>=0)
			{
				k = den-1;
				while(k>=1 && numExp[j]>=1)
				{
					if(denomExp[k]>=1 && (j+1)%(k+1)==0)
					{
						denomExp[k]--;
						numExp[j]--;
						numExp[((j+1)/(k+1))-1]++;
					}
					else k--;
				}
				j--;
			}
		}
		else if(numExp.length==0)
		{
			denomExp = new int[1];
			denomExp[0] = 1;
		}
	}
	
	public Fraction multiply(Fraction f)
	{
		int[] num = new int[0];
		int[] den = {1};
		if(this.numExp.length>=0 && f.numExp.length>=0 && this.denomExp.length>=0 && f.denomExp.length>=0)
		{
			num = new int[Math.max(this.numExp.length,f.numExp.length)];
			den = new int[Math.max(this.denomExp.length,f.denomExp.length)];
			for(int i=0;i<Math.max(num.length, den.length);i++)
			{
				if(i<this.numExp.length) num[i]+=this.numExp[i];
				if(i<f.numExp.length) num[i]+=f.numExp[i];
				if(i<this.denomExp.length) den[i]+=this.denomExp[i];
				if(i<f.denomExp.length) den[i]+=f.denomExp[i];
			}
		}		
		return new Fraction(num,den);
	}
	
	
	public static Fraction binomial(int n, int k)
	{
		Fraction binomial = null;
		if(n==0&&k==0)
			binomial = new Fraction(1,1);
		else if(n>=k&&n>=0&&k>=0)
		{
			int[] num = new int[n];
			int[] den = new int[Math.max(k, n-k)];
			for(int i=0;i<n;i++)
			{
				num[i]++;
				if(i<k) den[i]++;
				if(i<n-k) den[i]++;
			}
			binomial = new Fraction(num,den);
		}
		else
		{
			binomial = new Fraction(0,1);
		}
		return binomial;
	}
	
	public static int binValue(int n, int k)
	{
		return binomial(n,k).value();
	}
	
	public static Fraction simpleFrac(int n, int d)
	{
		return new Fraction(n,d);
	}
	
	public static void main(String[] args)
	{
		/*
		Fraction a = binomial(4,2);
		Fraction b = binomial(10,2);
		Fraction c = new Fraction(1,10);
		System.out.println(a.value());
		System.out.println(b.value());
		System.out.println(c.value());
		System.out.println(a.multiply(b).multiply(c).value());
		*/
		 //Fraction test = binomial(1,1);
		for(int i=-1;i<11;i++)
		{
			for(int j=-1;j<11;j++)
			{
				System.out.print(binomial(i,j).value()+",");
			}
			System.out.println();
		}
		//System.out.println(binomial(-1,-1).value());
		 //System.out.println(test.value());
		
	}
}
