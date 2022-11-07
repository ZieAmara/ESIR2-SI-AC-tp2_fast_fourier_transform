import java.io.IOException;

public class FFT_2D {

	//renvoie la TFD d'une image de complexes
	public static CpxImg FFT(CpxImg I) {

		CpxImg out = new CpxImg(I.taille());

		// FFT 1D sur les lignes
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(I.get_line(k)));
		  
		// transposition
		out.transpose();

		// FFT 1D sur les "nouvelles" lignes de out (anciennes colonnes)
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(out.get_line(k)));

		//on re transpose pour revenir dans le sens de départ
		out.transpose();
		
		//on divise par la taille de I
		out.multiply(1./I.taille());
		return out.recentrage();
	}
	
	//renvoie la TFD inverse d'une images de complexes
	public static CpxImg FFT_inverse(CpxImg I) {
		I = I.recentrage();
		CpxImg out = new CpxImg(I.taille());
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, I.get_line(k).conjugue());

		out = FFT(out).recentrage();
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, out.get_line(k).conjugue());
		return out;
	}

	// compression par mise é zéro des coefficients de fréquence 
	// FI contient la TDF de I 
	// Dans FI on met é zéros tous les coefficients correspondant é des fréquences inférieures é k
	public static void compression(CpxImg FI, int k) {
		// A COMPLETER
	}

	// compression par seuillage des coefficients faibles
	// FI contient la TDF de I 
	// Dans FI on met é zéros tous les coefficients dont le module est inférieur é seuil 
	// on renvoie le nombre de coefficients conservés 
	public static int compression_seuil(CpxImg FI, double seuil){
		//A COMPLETER
		return 0;
	}

	
	public static void main(String[] args) {
		
		try {			
			//PLACEZ ICI VOS TESTS en 2D
			//Exemple, lecture
			//BytePixmap BP = new BytePixmap("nomfichier.pgm");
			//CpxImg I = new CpxImg(BP);
			
			//Exemple, écriture
			//BP = I.convert_to_BytePixmap();
			//BP.write("nomfichier2.pgm");

				/*exo2 */
				BytePixmap bp_tigre_512 = new BytePixmap("tigre_512.pgm");	
				CpxImg I_tigre_512 = new CpxImg(bp_tigre_512);
				CpxImg fft_tigre_512 = FFT(I_tigre_512);
				CpxImg fft_inv_tigre_512 = FFT_inverse(fft_tigre_512);
				BytePixmap bp_tigre = fft_inv_tigre_512.convert_to_BytePixmap();
				bp_tigre.write("tigre.pgm");
				//Pixmap("tigre.pgm");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
