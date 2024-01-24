
TH1D* ReadROOTFileAsTh1D(const char* files[], 
	const int nFiles, 
	const char* histo_name,
	int bin_number,
	double min_value,
	double max_value)
{
	// histograms
	TH1D* h1 = new TH1D(Form("h1_%d", 1), histo_name, bin_number, min_value, max_value);

	for (int i = 0; i < nFiles; ++i)
	{
		TFile* file = TFile::Open(files[i]);
		TTree* tree = (TTree*)file->Get("variabletree");

		float_t energy = 0; //1.

		tree->SetBranchAddress("energy", &energy); //2.
		Long64_t n = tree->GetEntries();

		

		for (Long64_t i = 0; i < n; i++) {
			tree->GetEntry(i);
			h1->Fill(energy);
		}
		cout << n << endl;
		tree->ResetBranchAddresses();
	}
	TCanvas* c = new TCanvas();
	h1->Draw();
	return h1;
}

TH2D* ReadROOTFileAsTh2D(const char* files[],
	const int nFiles,
	const char* histo_name,
	int bin_number1,
	double min_value1,
	double max_value1,
	int bin_number2,
	double min_value2,
	double max_value2)
{
	// histograms
	TH2D* h1 = new TH2D(Form("h1_%d", 1), histo_name, bin_number1, min_value1, max_value1, bin_number2, min_value2, max_value2);

	for (int i = 0; i < nFiles; ++i)
	{
		TFile* file = TFile::Open(files[i]);
		TTree* tree = (TTree*)file->Get("variabletree");

		float_t energy = 0; //1.
		float_t MC_energy = 0;

		tree->SetBranchAddress("energy", &energy); //2.
		tree->SetBranchAddress("kinE_start_MC", &MC_energy);
		Long64_t n = tree->GetEntries();

		for (Long64_t i = 0; i < n; i++) {
			tree->GetEntry(i);
			h1->Fill(energy, MC_energy);
		}
		cout << n << endl;
		tree->ResetBranchAddresses();
	}
	TCanvas* c = new TCanvas();
	h1->Draw("COLZ");
	return h1;
}

void StretchHistogram(TH1D* originalHist, double scaleFactor) {
	// 获取原始直方图的参数
	int nBins = originalHist->GetNbinsX();
	double oldXmin = originalHist->GetXaxis()->GetXmin();
	double oldXmax = originalHist->GetXaxis()->GetXmax();

	// 计算新的轴范围
	double newXmin = oldXmin * scaleFactor;
	double newXmax = oldXmax * scaleFactor;

	// 创建新的直方图
	TH1D* stretchedHist = new TH1D("stretchedHist", originalHist->GetTitle(),
		nBins, newXmin, newXmax);

	// 复制 bin 的内容
	for (int i = 1; i <= nBins; ++i) {
		stretchedHist->SetBinContent(i, originalHist->GetBinContent(i));
		stretchedHist->SetBinError(i, originalHist->GetBinError(i));
	}

	// 可以选择将新的直方图返回或进行其他操作
	// ...
}


void unfold()
{
	const int binNum = 8;

	// measured spectrum
	const char* files_exp[] = { "./inputFiles/20231222132245574_bkg_1h_Vfc829_Vm861_Va1251_Vm_anti98_Va_anti488_25MHz_500ns_1pC_tree.root",
	"./inputFiles/20231222142844335_bkg_1h_Vfc829_Vm861_Va1251_Vm_anti98_Va_anti488_25MHz_500ns_1pC_tree.root" };
	const int nFiles_exp = 2;
	const char* histo_name_exp = "Measured energy - experiment; Measured Channel; Counts";
	TH1D* h_measured_energy = ReadROOTFileAsTh1D(files_exp,
		nFiles_exp,
		histo_name_exp, binNum, 0, 50000);
	//h_measured_energy->Rebin(6);
	//StretchHistogram(h_measured_energy, 0.9);

	// observed spectrum
	const char* files_obs[] = { "./inputFiles/Rawroot_gamma_digi_0_tree.root",
	"./inputFiles/Rawroot_gamma_digi_1_tree.root" };
	const int nFiles_obs = 2;
	const char* histo_name_obs = "Observed energy - simulation; Observed Channel; Counts";
	TH1D* h_observed_energy = ReadROOTFileAsTh1D(files_obs,
		nFiles_obs,
		histo_name_obs, binNum, 0, 50000);

	// initial spectrum
	TH1D* h_initial_energy = new TH1D("h1", "Initial energy - simulation;initial energy;counts", 256, 0, 3.044);
	std::ifstream file("./inputFiles/unfoldedSpectra_centerValue.txt");
	std::string line;
	int i = 1;
	while (std::getline(file, line)) 
	{
		double value = std::stod(line);
		h_initial_energy->SetBinContent(i, value);
		i++;
	}
	file.close();
	double currentTotal = h_initial_energy->Integral();
	double scaleFactor = 8e8 / currentTotal;
	h_initial_energy->Scale(scaleFactor);
	h_initial_energy->Rebin(256 / binNum);
	TCanvas* c = new TCanvas();
	h_initial_energy->Draw();

	// response matrix
	const char* histo_name_response = "response; Observed energy; True energy";
	TH2D* h_response = ReadROOTFileAsTh2D(files_obs,
		nFiles_obs,
		histo_name_response, binNum, 0, 50000, binNum, 0, 3.044);
	//h_response->Rebin2D(6, 6);

	// unfolding
	TH1D* nullHist = new TH1D("nullHist", "", binNum, 0, 3.044);
	RooUnfoldResponse response(nullHist, h_initial_energy, h_response, 0, 0);
	RooUnfoldBayes unfold_bayes(&response, h_measured_energy, 1000);
	TH1* hReco_bayes = unfold_bayes.Hunfold();
	TH1D* h_unfolded = new TH1D("", "h_unfolded;Energy;Counts", binNum, 0, 3.044);
	for (int i = 1; i <= 6; i++)
	{
		cout << Form("%d-th bin:\n",i);
		cout << hReco_bayes->GetBinContent(i) << endl;
		cout << hReco_bayes->GetBinError(i) << endl;
		h_unfolded->SetBinContent(i, hReco_bayes->GetBinContent(i));
	}

	TCanvas* c_result = new TCanvas();
	c_result->cd();
	h_unfolded->Draw();
}