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
				BytePixmap bp_tigre_512 = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");	
				CpxImg I_tigre_512 = new CpxImg(bp_tigre_512);
				CpxImg fft_tigre_512 = FFT(I_tigre_512);
				CpxImg fft_inv_tigre_512 = FFT_inverse(fft_tigre_512);
				BytePixmap bp_tigre = fft_inv_tigre_512.convert_to_BytePixmap();
				bp_tigre.write("AC_tp2_part2_donnees/tigre.pgm");

				/*exo3 */
				//mire1
				BytePixmap bp_mire1 = new BytePixmap("AC_tp2_part2_donnees/mire1.pgm");
				CpxImg fft_mire1 = FFT(new CpxImg(bp_mire1));
				bp_mire1 = fft_mire1.convert_to_BytePixmap();
				bp_mire1.write("AC_tp2_part2_donnees/mire11.pgm");
				//mire2
				BytePixmap bp_mire2 = new BytePixmap("AC_tp2_part2_donnees/mire2.pgm");
				CpxImg fft_mire2 = FFT(new CpxImg(bp_mire2));
				bp_mire2 = fft_mire2.convert_to_BytePixmap();
				bp_mire2.write("AC_tp2_part2_donnees/mire22.pgm");
				//mire3
				BytePixmap bp_mire3 = new BytePixmap("AC_tp2_part2_donnees/mire3.pgm");
				CpxImg fft_mire3 = FFT(new CpxImg(bp_mire3));
				bp_mire3 = fft_mire3.convert_to_BytePixmap();
				bp_mire3.write("AC_tp2_part2_donnees/mire33.pgm");
				//fingerprint
				BytePixmap bp_fingerprint = new BytePixmap("AC_tp2_part2_donnees/fingerprint.pgm");
				CpxImg fft_fingerprint = FFT(new CpxImg(bp_fingerprint));
				bp_fingerprint = fft_fingerprint.convert_to_BytePixmap();
				bp_fingerprint.write("AC_tp2_part2_donnees/fingerprint1.pgm");
				//barbara
				BytePixmap bp_barbara_512 = new BytePixmap("AC_tp2_part2_donnees/barbara_512.pgm");
				CpxImg fft_barbara_512 = FFT(new CpxImg(bp_barbara_512));
				bp_barbara_512 = fft_barbara_512.convert_to_BytePixmap();
				bp_barbara_512.write("AC_tp2_part2_donnees/barbara1.pgm");

				/*exo4 */
				//
				

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
