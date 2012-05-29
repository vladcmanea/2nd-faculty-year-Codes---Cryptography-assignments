#include "RSA.h"
#include "CPrecisionTimer.h"

int main()
{
	try
	{
		// create the system pointers
		RSA *cs;
		CPrecisionTimer *timer;
		double elapsed;
		bool x;

		// Wiener or not?
		cout << "Write 1 to cope with Wiener and 0 not to do so\n";
		cin >> x;

		// prompt for errors
		while (x != 0 && x != 1)
		{
			cout << "\nWrite 1 to cope with Wiener and 0 not to do so\n";
			cin >> x;
		}

		timer = new CPrecisionTimer();
		
		timer->Start();
		cout << "\nRSA parameters generating...\n";

		// generate output format
		ofstream fp1("public.txt");
		ofstream fs1("secret.txt");
		cs = new RSA(fp1, fs1, x);
		fp1.close();
		fs1.close();

		elapsed = timer->Stop();
		cout << "RSA parameters generated in " << elapsed << " seconds\n";

		timer->Start();
		cout << "\nRSA plain encrypting...\n";

		// encrypt content
		ifstream fc2("plain.txt");
		ifstream fp2("public.txt");
		ofstream fa2("encrypted.txt");
		cs->applyEncrypt(fc2, fp2, fa2);
		fc2.close();
		fp2.close();
		fa2.close();

		elapsed = timer->Stop();
		cout << "RSA plain encrypted in " << elapsed << " seconds\n";

		timer->Start();
		cout << "\nRSA encrypted decrypting...\n";

		// decrypt content
		ifstream fc3("encrypted.txt");
		ifstream fp3("public.txt");
		ifstream fs3("secret.txt");
		ofstream fa3("decrypted.txt");
		cs->applyDecrypt(fc3, fp3, fs3, fa3);
		fc3.close();
		fp3.close();
		fs3.close();
		fa3.close();

		elapsed = timer->Stop();
		cout << "RSA encrypted decrypted in " << elapsed << " seconds\n";

		// just in Wiener case
		if (x == 0)
		{
			timer->Start();
			cout << "\nRSA encrypted decrypting using Wiener...\n";

			// decrypt content by using Wiener's attack
			ifstream fc4("encrypted.txt");
			ifstream fp4("public.txt");
			ofstream fa4("decryptedWiener.txt");
			cs->applyWienerDecrypt(fc4, fp4, fa4);
			fc4.close();
			fp4.close();
			fa4.close();

			elapsed = timer->Stop();
			cout << "RSA encrypted decrypted using Wiener in " << elapsed << " seconds\n";
		}

		timer->Start();
		cout << "\nRSA factoring...\n";

		// factor content
		ifstream fp5("public.txt");
		ifstream fs5("secret.txt");
		ofstream fa5("factors.txt");
		cs->applyFactors(fp5, fs5, fa5);
		fp5.close();
		fs5.close();
		fa5.close();

		elapsed = timer->Stop();
		cout << "RSA factored in " << elapsed << " seconds\n";

		cout << '\n';

		// delete the system pointers
		delete cs;
		delete timer;
	}
	catch (char* e)
	{
		// write error
		cout << e << '\n';
	}

	return 0;
}