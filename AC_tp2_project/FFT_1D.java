import java.lang.Math;

/**** on va ici implémenter la transformée de Fourier rapide 1D ****/

public class FFT_1D {
	
	// "combine" c1 et c2 selon la formule vue en TD
	// c1 et c2 sont de même taille
	// la taille du résultat est le double de la taille de c1
	public static CpxTab combine(CpxTab c1, CpxTab c2) {
		assert (c1.taille()==c2.taille()) : "combine: c1 et c2 ne sont pas de même taille, taille c1="+c1.taille()+" taille c2="+c2.taille();
		//A FAIRE
    int n=2*c1.taille(); // On stock la taille de c1 
    /* CpxTab c va stocker la combinaison de c1 et c2
     * 		c[k] 		 = c1[k] + w[k]*c2[k]
		 * 		c[k+n/2] = c1[k] - w[k]*c2[k]
		 * Ces formules ont servi pour le calcul des parties réelles et imaginaire de c
		 * après calcul on obtient :
		 * 		Reel(c[k]) 		 = Reel(c1[k]) + cos((2kPI)/n)*Reel(c2[k]) - sin((2kPI)/n)*Im(c2[k])
		 * 		Im(c[k])			 = Im(c1[k]) 	 + sin((2kPI)/n)*Reel(c2[k]) + cos((2kPI)/n)*Im(c2[k])
		 * 		Reel(c[k+n/2]) = Reel(c1[k]) - cos((2kPI)/n)*Reel(c2[k]) + sin((2kPI)/n)*Im(c2[k])
		 * 		Im(c[k+n/2])	 = Im(c1[k]) 	 - sin((2kPI)/n)*Reel(c2[k]) - cos((2kPI)/n)*Im(c2[k])
     * */
    CpxTab c = new CpxTab(n); 
    // La boucle for suivante permet de remplir c
    for (int k=0; k<n/2; k++){
      // on calcul les partie réelles et imaginaire de c
      double c_k_reel = c1.get_p_reel(k)+Math.cos((2*k*Math.PI)/n)*c2.get_p_reel(k)-Math.sin((2*k*Math.PI)/n)*c2.get_p_imag(k);
      double c_k_imag = c1.get_p_imag(k)+Math.sin((2*k*Math.PI)/n)*c2.get_p_reel(k)+Math.cos((2*k*Math.PI)/n)*c2.get_p_imag(k);
      double c_kn2_reel = c1.get_p_reel(k)-Math.cos((2*k*Math.PI)/n)*c2.get_p_reel(k)+Math.sin((2*k*Math.PI)/n)*c2.get_p_imag(k);
      double c_kn2_imag = c1.get_p_imag(k)-Math.cos((2*k*Math.PI)/n)*c2.get_p_imag(k)-Math.sin((2*k*Math.PI)/n)*c2.get_p_reel(k);
      // on stosk ces partie réelles et imaginaire dans c
      c.set_p_reel(k, c_k_reel);
      c.set_p_imag(k, c_k_imag);
      c.set_p_reel(k+n/2, c_kn2_reel);
      c.set_p_imag(k+n/2, c_kn2_imag);
    }
		return c;
	}

	//renvoie la TFD d'un tableau de complexes
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(CpxTab x) {
		//A FAIRE : Test d'arrêt
    if (x.taille()==1) return x;

		assert (x.taille()%2==0) : "FFT: la taille de x doit être une puissance de 2";
		
		//A FAIRE : Décomposition en "pair" et "impair" et appel récursif
    int n = x.taille()/2;
		CpxTab x_pair = new CpxTab(n);
    CpxTab x_impair = new CpxTab(n);
    int i=0;
		int j=0;
		while (i<n || j<n) {
			for (int k=0; k<2*n; k++){
				if (k%2==0){
					x_pair.set_p_reel(i, x.get_p_reel(k));
					x_pair.set_p_imag(i, x.get_p_imag(k));
					i++;
				}
				else {
					x_impair.set_p_reel(j, x.get_p_reel(k));
					x_impair.set_p_imag(j, x.get_p_imag(k));
					j++;
				}
			}
		}
		CpxTab y_pair = FFT(x_pair);
		CpxTab y_impair = FFT(x_impair);
		CpxTab y = combine(y_pair, y_impair);
		return y;
	}

	//renvoie la TFD d'un tableau de réels
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(double[] x) {
		return FFT(new CpxTab(x));
	}
			
	//renvoie la transformée de Fourier inverse de y
	public static CpxTab FFT_inverse(CpxTab y) {
		//A FAIRE
		int n = y.taille();
		CpxTab c1 = (FFT(y.conjugue())).conjugue();
		CpxTab y_FFT_inverse = new CpxTab(n);
		for (int i=0; i<n; i++){
			y_FFT_inverse.set_p_reel(i, c1.get_p_reel(i)/n);
			y_FFT_inverse.set_p_imag(i, c1.get_p_imag(i)/n);
		}
		return y_FFT_inverse;
	}
	
	//calcule le produit de deux polynômes en utilisant la FFT
	//tab1 et tab2, sont les coefficients de ces polynômes
	// CpxTab sera le tableau des coefficients du polynôme produit (purement réel)
	public static CpxTab multiplication_polynome_viaFFT(double[] tab1, double[] tab2) {
		
		//on commence par doubler la taille des vecteurs en rajoutant des zéros à la fin (cf TD)
		double[] t1 = new double[2*tab1.length], t2 = new double[2*tab1.length];
		for (int i = 0; i < tab1.length; i++) {
			t1[i] = tab1[i];
			t2[i] = tab2[i];
		}
		CpxTab fft_tab1=FFT(t1);
		CpxTab fft_tab2=FFT(t2);
		//On fait la multiplication par point des deux polynomes representés par point
		CpxTab produit = CpxTab.multiplie(fft_tab1, fft_tab2);
		CpxTab fft_inv_produit = FFT_inverse(produit);
		//On retrouve ensuite les coefficients du polynomes produits en faisant la FFT inverse du résultat obtenu
		return fft_inv_produit;
	}

	
	//renvoie un tableau de réels aléatoires
	//utile pour tester la multiplication de polynômes
	public static double[] random(int n) {
		double[] t = new double[n];

		for (int i = 0; i < n; i++)
			t[i] = Math.random();
		return t;
	}

	//effectue la multiplication de polynômes représentés par coefficients
	// p1, p2 les coefficients des deux polynômes P1 et P2
	// renvoie les coefficients du polynôme P1*P2
	private static double [] multiplication_polynome_viaCoeff(double[] p1, double[] p2){
		
		int n = p1.length + p2.length - 1;
		double a,b;
		double [] out = new double[n];
		for (int k = 0; k < n; k++) {
			for (int i = 0; i <= k; i++) {
				a = (i<p1.length) ? p1[i]:0;
				b = (k-i<p2.length) ? p2[k-i] : 0;
				out[k] += a*b;
			}
		}
		return out;
	}
	

	//affiche un tableau de réels
	private static void afficher(double [] t){
		System.out.print("[");
		for(int k=0;k<t.length;k++){
			System.out.print(t[k]);
			if (k<(t.length-1))
				System.out.print(" ");
		}
		System.out.println("]");
	}
	
	public static void main(String[] args) {
		double[] t5 = {1,2,3,4};

		/* Exo 2: calculez et affichez TFD(1,2,3,4) */
			//A FAIRE
			System.out.println("-----------------------------------------------------");
			System.out.println("calcule et affichage de TFD(1,2,3,4)");
			System.out.println("	" + FFT(t5).toString());
		
		/* Exo 3: calculez et affichez TFD_inverse(TFD(1,2,3,4)) */
			//A FAIRE		
			System.out.println("-----------------------------------------------------");
			System.out.println("calcule et affichage de TFD_inverse(TFD(1,2,3,4))");
			System.out.println("	" + FFT_inverse(FFT(t5)).toString());

		/* Exo 4: multiplication polynomiale, vérification*/
			/* A(X) = 2 et B(X)=-3 */
			//A FAIRE		
			System.out.println("-----------------------------------------------------");
			System.out.println("Comparaison des 2 méthodes de multiplications polynomiales\n\tPour A(X) = 2 et B(X)=-3 :");
			double[] t7 = {2};
			double[] t8 = {-3};
			System.out.println("\t\tmult via FFT  --> " + multiplication_polynome_viaFFT(t7, t8));
			System.out.print(  "\t\tmult via coeff -> ");
			afficher(multiplication_polynome_viaCoeff(t7, t8));

			/* A(X) = 2+X et B(X)= -3+2X */
			//A FAIRE
			System.out.println("\tPour A(X) = 2+X et B(X)= -3+2X :");
			double[] t9 = {2, 1};
			double[] t10 = {-3,2};
			System.out.println("\t\tmult via FFT  --> " + multiplication_polynome_viaFFT(t9, t10));
			System.out.print(  "\t\tmult via coeff -> ");
			afficher(multiplication_polynome_viaCoeff(t9, t10));

			/* A(X) = 1 + 2X + 3X^2 + 4X^3 et B(X) = -3 + 2X - 5 X^2*/
			System.out.println("\tPour A(X) = 1 + 2X + 3X² + 4X³ et B(X) = -3 + 2X - 5 X² :");
			double[] t6 = {-3,2,-5,0};
			System.out.println("\t\tmult via FFT  --> " + multiplication_polynome_viaFFT(t5, t6));
			System.out.print(  "\t\tmult via coeff -> ");
			afficher(multiplication_polynome_viaCoeff(t5, t6));
	

		/* Exo 5: comparaison des temps de calculs */
	
		// Pour étude du temps de calcul 
		int[] t = {124, 250, 500, 1000, 2000};  // taille des polynômes à multiplier (testez différentes valeurs en gardant des puissances de 2)
		System.out.println("-----------------------------------------------------");
		for (int n: t){
			System.out.println("Temps de calcul pour n="+n);
			double[] tab1 = random(n), tab2 = random(n);
			long date1, date2;
			date1 = System.currentTimeMillis();
			multiplication_polynome_viaCoeff(tab1, tab2);
			date2 = System.currentTimeMillis();
			System.out.println("   via Coeff: " + (date2 - date1));

			date1 = System.currentTimeMillis();
			multiplication_polynome_viaFFT(tab1, tab2);
			date2 = System.currentTimeMillis();
			System.out.println("   via FFT  : " + (date2 - date1));
		}
	
	}

}