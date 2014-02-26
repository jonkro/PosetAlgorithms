import java.util.HashSet;


public class Element 
{
	HashSet<Integer> upShadow;
	HashSet<Integer> downShadow;
	int id;
	
	Element(int i)
	{
		upShadow = new HashSet<Integer>();
		downShadow = new HashSet<Integer>();
		id = i;
	}
	
	Element(int i, HashSet<Integer> u, HashSet<Integer> d)
	{
		upShadow = u;
		downShadow = d;
		id = i;
	}
	
	public Element hardCopy(int ident)
	{
		Element copy = new Element(ident);
		for(int i : upShadow)
		{
			copy.newGreater(i);
		}
		for(int i : downShadow)
		{
			copy.newSmaller(i);
		}
		return copy;
	}
	
	public void deleteGreater(Integer g)
	{
		upShadow.remove(g);
	}
	
	public void deleteSmaller(Integer s)
	{
		downShadow.remove(s);
	}
	

	
	void newSmaller(Integer i)
	{
		downShadow.add(i);
	}
	
	void newGreater(Integer i)
	{
		upShadow.add(i);
	}
	
	int getID()
	{
		return id;
	}
	
	HashSet<Integer> getUpShadow()
	{
		return upShadow;
	}
	
	HashSet<Integer> getDownShadow()
	{
		return downShadow;
	}
	
	boolean isMinumum()
	{
		//System.out.println(downShadow.isEmpty());
		if(downShadow.isEmpty()) return true;
		else return false;
	}
	
	boolean isMaximum()
	{
		if(upShadow.isEmpty()) return true;
		else return false;
	}
	
	public Element copy()
	{
		HashSet<Integer> u = new HashSet<Integer>(); 
		HashSet<Integer> d = new HashSet<Integer>(); 
		for(int i : upShadow)
		{
			u.add(i);
		}
		for(int i : downShadow)
		{
			d.add(i);
		}
		return new Element(this.getID(),u,d);
	}
	
	public Element copy(int offset)
	{
		HashSet<Integer> u = new HashSet<Integer>(); 
		HashSet<Integer> d = new HashSet<Integer>(); 
		for(int i : upShadow)
		{
			u.add(i+offset);
		}
		for(int i : downShadow)
		{
			d.add(i+offset);
		}
		return new Element(this.getID()+offset,u,d);
	}
	
	public String toString()
	{
		StringBuffer str = new StringBuffer();
		str.append("[");
		for(Integer i : downShadow)
		{
			str.append(i+",");
		}
		str.append("]->"+this.getID()+"->[");
		for(Integer i : upShadow)
		{
			str.append(i+",");
		}
		str.append("]");
		return str.toString();
	}
}