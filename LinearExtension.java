import java.util.Vector;


public class LinearExtension 
{
	Vector<Integer> linExt;
	int jumps;
	
	public LinearExtension()
	{
		linExt = new Vector<Integer>();
		jumps = 0;
	}
	
	public LinearExtension(LinearExtension original)
	{
		linExt = new Vector<Integer>();
		for(int i=0; i<original.size(); i++)
		{
			linExt.add(original.getElemAt(i));
		}
		jumps = original.getJumps();
	}
	
	public void addJump()
	{
		jumps++;
	}
	
	public void addElem(int e)
	{
		linExt.add(e);
	}
	
	public int getJumps()
	{
		return jumps;
	}
	
	public int size()
	{
		return linExt.size();
	}
	
	public int getElemAt(int i)
	{
		return linExt.get(i);
	}
	
	public int lastElem()
	{
		return linExt.lastElement().intValue();
	}
	
	public void add(int position, int elem)
	{
		linExt.add(position, elem);
	}
	
	public void add(int elem)
	{
		linExt.add(elem);
	}
	
	public void add(int position, int elem, boolean jump)
	{
		add(position,elem);
		if(jump) addJump();
	}
	
	public void add(int elem, boolean jump)
	{
		add(elem);
		if(jump) addJump();
	}

	public boolean isEmpty() 
	{
		return linExt.isEmpty();
	}
	
	public boolean jumpBefore(int i)
	{
		return !((linExt.get(i)-linExt.get(i-1)==1)||(linExt.get(i)-linExt.get(i-1)==(linExt.size()/2)));
	}
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		for(int i=0; i<linExt.size()-1;i++)
		{
			str.append(linExt.get(i));
			str.append(",");
		}
		str.append(linExt.get(linExt.size()-1));
		return str.toString();
	}
	
	public String showMiddleChain()
	{
		StringBuffer str = new StringBuffer();
		int j = 0;
		for(int i=1; i<linExt.size()-1;i++)
		{
			if(jumpBefore(i)) j++;
			if(j==2) 
			{
				str.append("(");
			}
			str.append(linExt.get(i));
			if(j==2) 
			{
				str.append(")");
			}

				str.append(",");
		}
		return str.toString();
	}
}
