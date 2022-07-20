void cutData() {

    TFile *output = new TFile("muons_peak.root", "RECREATE");
    TTree *tree = new TTree("tree", "tree");
    
    double t, mass_inv;
    int id, count, r;
    
    tree->Branch("time", &t, "t/D");
    tree->Branch("minv", &mass_inv, "minv/D");
    tree->Branch("id",&id,"id/I");
    tree->Branch("count",&count,"count/I");
    tree->Branch("run",&r,"r/I");
    tree->SetAutoSave(150000000);
    
    
    TFile *f = new TFile("data.root");
    TTreeReader reader("events",f);
    TTreeReaderValue<int>ev(reader,"event");
    TTreeReaderValue<double>time(reader,"time");
    TTreeReaderValue<int>run(reader,"run");
    TTreeReaderValue<Double_t>mu0_pid(reader,"mu0_pid");
    TTreeReaderValue<Double_t>mu1_pid(reader,"mu1_pid");
    TTreeReaderValue<TVector3>pMu0(reader,"mu0_p");
    TTreeReaderValue<TVector3>pMu1(reader,"mu1_p");
    
    
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,9);
    
    
    
    const double m_mu = 0.105658;
    const double m_mu_2 = m_mu*m_mu;
    
    
    int i = 0;
    count = 0;
    while (reader.Next()) {
        double E0;
        double E1;
        double E;
        TVector3 p;
        double size_p_2; 
        double Minv;   
        
        if ((*mu0_pid < 0.9) or (*mu1_pid < 0.9)) continue;
                
            
        p = *pMu0 + *pMu1; //vector sum of momenta
        size_p_2 = p.Mag2(); //squared magnitude of vector p 
                
        E0 = sqrt(pMu0->Mag2() + m_mu_2); 
        E1 = sqrt(pMu1->Mag2() + m_mu_2);
        E = E0 + E1;
                
        Minv = sqrt(E*E - size_p_2); // inv hmota
                
        if (((Minv > 10.2 and (Minv < 10.8)) and (*run < 4000))) {
            id = distribution(generator); 
            mass_inv = Minv;
            t = *time;
            r = *run;
            count++;
            tree->Fill();
        }   
    }
    cout << count << endl;    
    output->Write();
    
    output->Close();

    


}
