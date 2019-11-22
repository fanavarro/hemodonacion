package dis.um.es.hemodonacionML.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

public class Utils {
	public static BufferedReader readDataFile(File file) throws FileNotFoundException{
		BufferedReader inputReader = null;
		inputReader = new BufferedReader(new FileReader(file));
		
		return inputReader;
	}
}
