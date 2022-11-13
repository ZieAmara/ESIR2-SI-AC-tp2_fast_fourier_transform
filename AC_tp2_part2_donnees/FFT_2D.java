import java.io.IOException;

public class FFT_2D {

	//renvoie la TFD d'une image de complexes
	public static CpxImg FFT(CpxImg I) {

		CpxImg out = new CpxImg(I.taille());

		// FFT 1D sur les lignes
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D_part2.FFT(I.get_line(k)));
		  
		// transposition
		out.transpose();

		// FFT 1D sur les "nouvelles" lignes de out (anciennes colonnes)
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D_part2.FFT(out.get_line(k)));

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
		/*exo 5 */
		int milieu = FI.taille()/2;
		for(int i=0; i<FI.taille(); i++) {
			for(int j=0; j<FI.taille(); j++) {
				if ((i < milieu-k) || (i > milieu+k)) {
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
				else if((j < milieu-k) || (j > milieu+k)){
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
			}
		}		
	}

	

	// compression par seuillage des coefficients faibles
	// FI contient la TDF de I 
	// Dans FI on met � z�ros tous les coefficients dont le module est inf�rieur � seuil 
	// on renvoie le nombre de coefficients conserv�s 
	public static int compression_seuil(CpxImg FI, double seuil){
		/*exo 6 */
		int k =0;
		for(int l=0; l<FI.taille(); l++) {
			for(int m=0; m<FI.taille(); m++) {
				if (Math.sqrt(FI.get_p_reel(l, m)*FI.get_p_reel(l, m) + FI.get_p_imag(l, m)*FI.get_p_imag(l, m)) < seuil) {
					FI.set_p_reel(l, m, 0);
					FI.set_p_imag(l, m, 0);
				}
				else k++;		
			}
		}
		System.out.println(k);
		return k;
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

			/*exo 2 */
				BytePixmap bp_tigre_512 = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");	
				CpxImg I_tigre_512 = new CpxImg(bp_tigre_512);
				CpxImg fft_tigre_512 = FFT(I_tigre_512);
				CpxImg fft_inv_tigre_512 = FFT_inverse(fft_tigre_512);
				BytePixmap bp_tigre = fft_inv_tigre_512.convert_to_BytePixmap();
				bp_tigre.write("AC_tp2_part2_donnees/tigre.pgm");

			/*exo 3 */
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

			/*exo 4 */
				// basse fréquences
				CpxImg fft_tigre_b = fft_tigre_512;
				fft_tigre_b.set_p_imag(fft_tigre_512.taille()/2, fft_tigre_512.taille()/2, 0);
				fft_tigre_b.set_p_reel(fft_tigre_512.taille()/2, fft_tigre_512.taille()/2, 0);
				bp_tigre = FFT_inverse(fft_tigre_b).convert_to_BytePixmap();
				bp_tigre.write("AC_tp2_part2_donnees/tigre_basse_fréquence.pgm");
				// basse fréquences +1
				CpxImg fft_tigre_b1 = fft_tigre_512;
				fft_tigre_b1.set_p_imag(fft_tigre_512.taille()/2+1, fft_tigre_512.taille()/2, 0);
				fft_tigre_b1.set_p_reel(fft_tigre_512.taille()/2+1, fft_tigre_512.taille()/2, 0);
				bp_tigre = FFT_inverse(fft_tigre_b1).convert_to_BytePixmap();
				bp_tigre.write("AC_tp2_part2_donnees/tigre_basse1_fréquence.pgm");
				// haute fréquences
				CpxImg fft_tigre_h = fft_tigre_512;
				fft_tigre_h.set_p_imag(fft_tigre_512.taille()/2, fft_tigre_512.taille()-1, 0);
				fft_tigre_h.set_p_reel(fft_tigre_512.taille()/2, fft_tigre_512.taille()-1, 0);
				bp_tigre = FFT_inverse(fft_tigre_h).convert_to_BytePixmap();
				bp_tigre.write("AC_tp2_part2_donnees/tigre_haute_fréquence.pgm");

				/*exo 7 */
				// k = sqrt(0.9/4)*n
				CpxImg fft_barbara_compression = fft_barbara_512;
				int k = (int) (fft_barbara_compression.taille()*Math.sqrt(0.1/4));
				System.out.println(fft_barbara_compression.taille()*fft_barbara_compression.taille()-4*k*k);
				compression(fft_barbara_compression, k);
				bp_barbara_512 = FFT_inverse(fft_barbara_compression).convert_to_BytePixmap();
				bp_barbara_512.write("AC_tp2_part2_donnees/barbara_comprimé.pgm");
				// compression seuil
				CpxImg fft_barbara_compression_seuil = fft_barbara_512;
				double s = 29.9;//-9.915;
				compression_seuil(fft_barbara_compression_seuil, s);
				bp_barbara_512 = FFT_inverse(fft_barbara_compression_seuil).convert_to_BytePixmap();
				bp_barbara_512.write("AC_tp2_part2_donnees/barbara_comprimé_seuil.pgm");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
