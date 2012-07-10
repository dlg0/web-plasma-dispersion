#include <Wt/WApplication>
#include <Wt/WBreak>
#include <Wt/WContainerWidget>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WText>
#include <Wt/WDoubleValidator>
#include <Wt/WValidator>
#include <Wt/WColor>
#include <Wt/WComboBox>
#include <Wt/WGroupBox>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <sstream>
#include <Wt/WBorderLayout>
#include <Wt/WImage>
#include <Wt/Chart/WCartesianChart>
#include <Wt/WStandardItemModel>
#include <boost/lexical_cast.hpp>
#include <Wt/Chart/WAxis>
#include <stixVars.hpp>
#include <constants.hpp>

class IonSpec
{
		public:
		IonSpec(int, Wt::WGroupBox*);
		~IonSpec(void);
		void setVisible(void);
		void setHidden(void);

		private:

		public:
		Wt::WLineEdit *SpecAMU_LineEdit, *SpecZ_LineEdit;
		Wt::WBreak *SpecBreak;
		Wt::WText *SpecText;
		int SpecNo;
};

IonSpec::IonSpec(int SpecNoIn, Wt::WGroupBox *IonSpecies_Container )
{
	SpecNo = SpecNoIn;
	std::stringstream TmpStr;
	TmpStr << "Ion Species " << SpecNoIn << " (AMU, Z)";
	Wt::WString TmpStrWt(TmpStr.str());
	SpecText = new Wt::WText(TmpStrWt,IonSpecies_Container);
	SpecText->setHidden(1);
	SpecAMU_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit = new Wt::WLineEdit("1",IonSpecies_Container);
	SpecZ_LineEdit->setHidden(1);
	SpecBreak = new Wt::WBreak(IonSpecies_Container);
	SpecBreak->setHidden(1);
}
IonSpec::~IonSpec(void)
{
}
void IonSpec::setVisible(void)
{
	SpecText->setHidden(0);
	SpecAMU_LineEdit->setHidden(0);
	SpecZ_LineEdit->setHidden(0);
	SpecBreak->setHidden(0);
}
void IonSpec::setHidden(void)
{
	SpecText->setHidden(1);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit->setHidden(1);
	SpecBreak->setHidden(1);
}

class DispersionApp : public Wt::WApplication
{
		public:
		    DispersionApp(const Wt::WEnvironment& env);

		private:
			Wt::WPushButton *Update_PushButton;
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_;
			Wt::WLineEdit *Density_m3_LineEdit, *BField_LineEdit, *Freq_Hz_LineEdit;
			Wt::WDoubleValidator *DensityValidator, *BFieldValidator, *FreqValidator;
			Wt::WText *DensityValidText, *BFieldValidText, *FreqValidText;

			Wt::WContainerWidget *w;
			Wt::WBorderLayout *layout;

			std::vector<IonSpec> IonSpecies;

			Wt::WGroupBox *West_Container, *Center_Container, *ParameterType_Container, *IonSpecies_Container, 
					*PointParameterType_Container,*MagneticParameterType_Container;
			Wt::WButtonGroup *ParameterType_ButtonGroup;
			Wt::WComboBox *nIonSpecies_ComboBox;

			// Scan
			int nX;
			Wt::WText *nX_Text;
			Wt::WLineEdit *nX_LineEdit;

			// Magnetic field scan
			float r0, b0, rMin, rMax, kPar;
			Wt::WText *r0_Text, *b0_Text, *rMin_Text, *rMax_Text, *kPar_Text;
			Wt::WLineEdit *r0_LineEdit, *b0_LineEdit, *rMin_LineEdit, *rMax_LineEdit, 
					*kPar_LineEdit;

			std::vector <float> x,b;

			Wt::Chart::WCartesianChart *b_CartesianChart, *a_CartesianChart;
			Wt::WStandardItemModel *b_model, *a_model;
			Wt::Chart::WDataSeries *a_series, *b_series; 

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
			void UpdateIonSpecies();
			void RevealParameterSettings();
			void UpdateCalculation();
			void UpdateBField();

			void AddChart(std::string title, std::string xTitle, std::string yTitle, 
							Wt::Chart::WCartesianChart *_chart, Wt::WStandardItemModel *_model);

			int nSpecies;

};

DispersionApp::DispersionApp(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			w = new Wt::WContainerWidget(root());
			layout = new Wt::WBorderLayout();
			w->setLayout(layout,Wt::AlignTop | Wt::AlignJustify);

			//layout->addWidget(new Wt::WImage("naruto.jpg"),Wt::WBorderLayout::Center);

			enum CalculationType { Hot = 1, Cold = 2 };

			West_Container = new Wt::WGroupBox("West");
			West_Container->setId("West");
			West_Container->resize(250,400);
			layout->addWidget(West_Container,Wt::WBorderLayout::West);

			Center_Container = new Wt::WGroupBox("Center");
			Center_Container->setId("Center");
			layout->addWidget(Center_Container,Wt::WBorderLayout::Center);

			Wt::WButtonGroup *HotOrCold_ButtonGroup = new Wt::WButtonGroup(West_Container);
			Wt::WRadioButton *HotOrCold_Button;

			HotOrCold_Button = new Wt::WRadioButton("Hot", West_Container);
			new Wt::WBreak(West_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Hot);

			HotOrCold_Button = new Wt::WRadioButton("Cold", West_Container);
			new Wt::WBreak(West_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Cold);

			HotOrCold_ButtonGroup->setCheckedButton(HotOrCold_ButtonGroup->button(Cold));

			nX_Text = new Wt::WText("nX",West_Container);
			nX_LineEdit = new Wt::WLineEdit("100",West_Container);
			nX = boost::lexical_cast<int>(nX_LineEdit->text().narrow());

			enum ParameterType { Point = 1, MagneticField = 2, Numerical = 3 };

			ParameterType_Container = new Wt::WGroupBox("Parameter Scan Type");
			layout->addWidget(ParameterType_Container,Wt::WBorderLayout::East);
			ParameterType_ButtonGroup = new Wt::WButtonGroup(ParameterType_Container);

			Wt::WRadioButton *ParameterType_Button;

			ParameterType_Button = new Wt::WRadioButton("Point", ParameterType_Container);
			new Wt::WBreak(ParameterType_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,Point);

			ParameterType_Button = new Wt::WRadioButton("MagneticField", ParameterType_Container);
			new Wt::WBreak(ParameterType_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,MagneticField);

			ParameterType_Button = new Wt::WRadioButton("Numerical", ParameterType_Container);
			new Wt::WBreak(ParameterType_Container);
			ParameterType_ButtonGroup->addButton(ParameterType_Button,Numerical);

			ParameterType_ButtonGroup->setCheckedButton(ParameterType_ButtonGroup->button(MagneticField));
			ParameterType_ButtonGroup->checkedChanged().connect(
							boost::bind(&DispersionApp::RevealParameterSettings,this));

			PointParameterType_Container = new Wt::WGroupBox("Point Parameter Settings",root());

			new Wt::WBreak(West_Container);
			Wt::WText *DensityText = new Wt::WText(West_Container);
			DensityText->setText("Density [1/m^3]");
			DensityValidator = new Wt::WDoubleValidator(1e1,1e25,West_Container);
			DensityValidator->setMandatory(1);
			Density_m3_LineEdit = new Wt::WLineEdit("2e19",West_Container);
			Density_m3_LineEdit->setValidator(DensityValidator);
			DensityValidText = new Wt::WText(West_Container);
			Density_m3_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,Density_m3_LineEdit,DensityValidText));

			new Wt::WBreak(West_Container);

			Wt::WText *BFieldText = new Wt::WText(PointParameterType_Container);
			BFieldText->setText("BField [T]");
			BFieldValidator = new Wt::WDoubleValidator(1e-5,5e5,PointParameterType_Container);
			BFieldValidator->setMandatory(1);
			BField_LineEdit = new Wt::WLineEdit("5",PointParameterType_Container);
			BField_LineEdit->setValidator(BFieldValidator);
			BFieldValidText = new Wt::WText(PointParameterType_Container);
			BField_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,BField_LineEdit,BFieldValidText));
			//new Wt::WBreak(West_Container);

			Wt::WText *FreqText = new Wt::WText(West_Container);
			FreqText->setText("Freq [Hz]");
			FreqValidator = new Wt::WDoubleValidator(1e1,5e12,West_Container);
			FreqValidator->setMandatory(1);
			Freq_Hz_LineEdit = new Wt::WLineEdit("30e6",West_Container);
			Freq_Hz_LineEdit->setValidator(FreqValidator);
			FreqValidText = new Wt::WText(West_Container);
			Freq_Hz_LineEdit->changed().connect(
					boost::bind(&DispersionApp::editedLineEdit,this,Freq_Hz_LineEdit,FreqValidText));
			new Wt::WBreak(West_Container);

			PointParameterType_Container->setHidden(1);

			MagneticParameterType_Container = new Wt::WGroupBox("Magnetic Field Scan Settings",root());
			MagneticParameterType_Container->setHidden(0);

			kPar_Text = new Wt::WText("kPar",West_Container);
			kPar_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			new Wt::WBreak(West_Container);

			r0_Text = new Wt::WText("r0",West_Container);
			r0_LineEdit = new Wt::WLineEdit("2",West_Container);
			new Wt::WBreak(West_Container);

			b0_Text = new Wt::WText("b0",West_Container);
			b0_LineEdit = new Wt::WLineEdit("5",West_Container);
			new Wt::WBreak(West_Container);

			rMin_Text = new Wt::WText("rMin",West_Container);
			rMin_LineEdit = new Wt::WLineEdit("0.1",West_Container);
			new Wt::WBreak(West_Container);

			rMax_Text = new Wt::WText("rMax",West_Container);
			rMax_LineEdit = new Wt::WLineEdit("4",West_Container);

			Update_PushButton = new Wt::WPushButton("Update",West_Container);
			Update_PushButton->clicked().connect(this, &DispersionApp::UpdateCalculation);

			IonSpecies_Container = new Wt::WGroupBox("Ion Species",root());

			Wt::WText *nIonSpecText = new Wt::WText("N Ion Species",IonSpecies_Container);
			nIonSpecies_ComboBox = new Wt::WComboBox(IonSpecies_Container);
			nIonSpecies_ComboBox->addItem("1");
			nIonSpecies_ComboBox->addItem("2");
			nIonSpecies_ComboBox->addItem("3");
			nIonSpecies_ComboBox->addItem("4");
			nIonSpecies_ComboBox->addItem("5");
			nIonSpecies_ComboBox->addItem("6");
			nIonSpecies_ComboBox->setCurrentIndex(0);
			nIonSpecies_ComboBox->changed().connect(
				boost::bind(&DispersionApp::UpdateIonSpecies,this));

			new Wt::WBreak(IonSpecies_Container);

			nSpecies = 0;

			a_CartesianChart = new Wt::Chart::WCartesianChart(Center_Container); 
			b_CartesianChart = new Wt::Chart::WCartesianChart(Center_Container); 

			a_model = new Wt::WStandardItemModel(nX,2);
			b_model = new Wt::WStandardItemModel(nX,2);

			a_series = new Wt::Chart::WDataSeries(1,Wt::Chart::LineSeries);
			b_series = new Wt::Chart::WDataSeries(1,Wt::Chart::LineSeries);

			AddChart(std::string("a"),std::string("x [m]"),std::string("b [T]"), a_CartesianChart, a_model);
			AddChart(std::string("b"),std::string("x [m]"),std::string("b [T]"), b_CartesianChart, b_model);

			a_CartesianChart->addSeries(*a_series);
			b_CartesianChart->addSeries(*b_series);

			UpdateCalculation();
}

void DispersionApp::AddChart (std::string title, std::string xTitle, std::string yTitle, 
	 Wt::Chart::WCartesianChart *_Chart, Wt::WStandardItemModel *_model )
{

	//int nX = x.size();
	//_model->setHeaderData(0, Wt::WString("x (m)"));
	//_model->setHeaderData(1, Wt::WString("b (T)"));

	for (unsigned i = 0; i < nX; ++i) {
	  _model->setData(i, 0, i);
	  _model->setData(i, 1, i);
	}

	_Chart->setModel(_model);
	_Chart->setType(Wt::Chart::ScatterPlot);
	//_Chart->setLegendEnabled(true);
	_Chart->axis(Wt::Chart::XAxis).setLocation(Wt::Chart::ZeroValue);
	_Chart->resize(600,300);
	_Chart->setPlotAreaPadding(200,Wt::Left);
	_Chart->setPlotAreaPadding(80,Wt::Bottom);
	_Chart->setLegendLocation(Wt::Chart::LegendOutside,Wt::Top,Wt::AlignRight);
	_Chart->setTitle(Wt::WString(title));

	Wt::Chart::WAxis &yAxis = _Chart->axis(Wt::Chart::YAxis); 
	Wt::Chart::WAxis &xAxis = _Chart->axis(Wt::Chart::XAxis); 
	yAxis.setTitle(yTitle);
	xAxis.setTitle(xTitle);

	_Chart->setXSeriesColumn(0);
	//_Chart->addSeries(s);
}

void DispersionApp::RevealParameterSettings()
{
		if(ParameterType_ButtonGroup->checkedId()==1)
		{
			PointParameterType_Container->setHidden(0);
		}

		if(ParameterType_ButtonGroup->checkedId()==2)
		{
			PointParameterType_Container->setHidden(1);
		}

		if(ParameterType_ButtonGroup->checkedId()==3)
		{
			PointParameterType_Container->setHidden(1);
		}
}

void DispersionApp::UpdateIonSpecies()
{

		int nSpecTmp = nIonSpecies_ComboBox->currentIndex();
		std::cout << nSpecTmp << std::endl;

		for(int s=0; s<nSpecies; s++)
		{
				IonSpecies[s].setHidden();
		}

		IonSpecies.clear();
		nSpecies = 0;

		for(int s=0; s<=nSpecTmp; s++)
		{
				IonSpec mySpecParams(s+1,IonSpecies_Container);	
				IonSpecies.push_back(mySpecParams);
				IonSpecies[s].setVisible();
				nSpecies++;
		}

}

void DispersionApp::editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText)
{
		if(LineEditIn->validate() != Wt::WValidator::Valid)
		{
				LineEditIn->decorationStyle().setBackgroundColor(Wt::WColor("#F79A9A"));
				ValidText->setText("Not Valid Input");
		}
		else
		{
			LineEditIn->decorationStyle().setBackgroundColor(Wt::WColor("#66FF66"));
			ValidText->setText("Good to go");
		}
}

void DispersionApp::greet()
{
		    greeting_->setText("Hello there, " + nameEdit_->text());
}

void DispersionApp::UpdateBField ()
{
			nX = boost::lexical_cast<int>(nX_LineEdit->text().narrow());
			rMin = boost::lexical_cast<float>(rMin_LineEdit->text().narrow());
			rMax = boost::lexical_cast<float>(rMax_LineEdit->text().narrow());
			b0 = boost::lexical_cast<float>(b0_LineEdit->text().narrow());
			r0 = boost::lexical_cast<float>(r0_LineEdit->text().narrow());

			x.resize(nX);
			b.resize(nX);

			for(int i=0;i<nX;i++)
			{
					x[i] = i*(rMax-rMin)/nX+rMin;
					b[i] = b0*r0/x[i];
			}

			
}

void DispersionApp::UpdateCalculation()
{
			UpdateIonSpecies();
			UpdateBField();

			std::vector<double> tmp;

			for(int i=0;i<nX;i++)
			{

				float _bMag = b[i];
				std::vector<PlasmaSpecies> AllSpecies;
				double _ionDensity = 0.0;
				//ions
				for(int s=0;s<nSpecies;s++)
				{
						double _z;
						try {
							_z = boost::lexical_cast<double>(IonSpecies[s].SpecZ_LineEdit->text().narrow());
							std::cout << "Success: SpecZ_LineEdit read as " << _z << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecZ_LineEdit could not be read" << std::endl;
						}

						double _amu;
						try {
					   		_amu = boost::lexical_cast<double>(IonSpecies[s].SpecAMU_LineEdit->text().narrow());
							std::cout << "Success: SpecAMU_LineEdit read as " << _amu << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: SpecAMU_LineEdit could not be read" << std::endl;
						}

						double _n;
						try {
					   		_n = boost::lexical_cast<double>(Density_m3_LineEdit->text().narrow());
							std::cout << "Success: Density_m3_LineEdit read as " << _n << std::endl;
						} catch(boost::bad_lexical_cast &) {
							std::cout << "Error: Density_m3_LineEdit could not be read" << std::endl;
						}

						AllSpecies.push_back(PlasmaSpecies(_z,_amu,_n,_bMag));
						_ionDensity += _z * _n;
						std::cout << "_ionDensity: " << _ionDensity << std::endl;

						std::cout << "AllSpecies[" << s << "] wp = " << AllSpecies[s].wp << std::endl;
						std::cout << "AllSpecies[" << s << "] wc = " << AllSpecies[s].wc << std::endl;

				}
				//electrons
				AllSpecies.push_back(PlasmaSpecies(-1.0,_me_mi,_ionDensity,_bMag));

				std::cout << "AllSpecies[e] wp = " << AllSpecies[nSpecies].wp << std::endl;
				std::cout << "AllSpecies[e] wc = " << AllSpecies[nSpecies].wc << std::endl;

				double _freq = boost::lexical_cast<double>(Freq_Hz_LineEdit->text().narrow());
				double _omega = 2.0 * _pi * _freq;
				StixVars stix(_omega, AllSpecies);

				std::cout << "Stix.R = " << stix.R << std::endl;
				std::cout << "Stix.L = " << stix.L << std::endl;
				std::cout << "Stix.D = " << stix.D << std::endl;
				std::cout << "Stix.S = " << stix.S << std::endl;
				std::cout << "Stix.P = " << stix.P << std::endl;

				tmp.push_back(stix.R);

				AllSpecies.clear();
			}

			// Update plots

			std::cout << "Updating plots ..." << std::endl;

			int nRows = a_model->rowCount();
			if(nRows<nX) a_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) a_model->removeRows(0,nRows-nX);
			std::cout << "nX: "<<nX<<" nRows: "<<a_model->rowCount()<<std::endl;
			for (unsigned i = 0; i < nX; ++i) {
			  a_model->setData(i, 0, x[i]);
			  a_model->setData(i, 1, tmp[i]);
			}

			if(nRows<nX) b_model->insertRows(nRows-1,nX-nRows);
			if(nRows>nX) b_model->removeRows(0,nRows-nX);
			for (unsigned i = 0; i < nX; ++i) {
			  b_model->setData(i, 0, x[i]);
			  b_model->setData(i, 1, b[i]);
			}
}

Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
		    return new DispersionApp(env);
}

int main(int argc, char **argv)
{
		    return Wt::WRun(argc, argv, &createApplication);
}
