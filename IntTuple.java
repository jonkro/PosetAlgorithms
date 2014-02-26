
public class IntTuple 
{
	int lower;
	int higher;
	
	IntTuple(int a, int b)
	{
		if(a<=b)
		{
			lower = a;
			higher = b;
		}
		else
		{
			lower = b;
			higher = a;
		}
	}
	
	int getLower()
	{
		return lower;
	}
	
	int getHigher()
	{
		return higher;
	}
	
	  @Override
	  public int hashCode() 
	  {
	    final int prime = 31;
	    int result = 1;
	    result = prime * result + (lower + 1)* higher;
	    return result;
	  }
	/*
	@Override
	public boolean equals(IntTuple it)
	{
		boolean result = false;
		if(lower == it.getLower() && higher == it.getHigher()) result = true;
		return result;
	}*/
	  
		
	@Override
	public boolean equals(Object it)
	{
		boolean result = false;
		if(lower == ((IntTuple) it).getLower() && higher == ((IntTuple) it).getHigher()) result = true;
		return result;
	}
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		str.append("["+getLower()+","+getHigher()+"]");
		return str.toString();
	}
}
