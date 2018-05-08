//compute hamming distance value between two strings 
public class  hamming
{
   private String compOne;
   private String compTwo;

   public hamming(String one, String two)
   {
	   compTwo = two;
       compOne = one;
   }
   public int getHammingDistance(){
   		if (compOne.length() != compTwo.length())
   		{
   			return -1;
   		}
        int counter = 0;
        for (int i = 0; i < compOne.length(); i++)
        {
        	if (compOne.charAt(i) != compTwo.charAt(i)) counter++;
        }
        return counter;
    }    
}


