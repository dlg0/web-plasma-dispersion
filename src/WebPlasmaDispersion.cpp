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
	SpecAMU_LineEdit = new Wt::WLineEdit(IonSpecies_Container);
	SpecAMU_LineEdit->setHidden(1);
	SpecZ_LineEdit = new Wt::WLineEdit(IonSpecies_Container);
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

class HelloApplication : public Wt::WApplication
{
		public:
		    HelloApplication(const Wt::WEnvironment& env);

		private:
			Wt::WLineEdit *nameEdit_;
			Wt::WText *greeting_;
			Wt::WLineEdit *Density_m3_LineEdit, *BField_LineEdit, *Freq_Hz_LineEdit;
			Wt::WDoubleValidator *DensityValidator, *BFieldValidator, *FreqValidator;
			Wt::WText *DensityValidText, *BFieldValidText, *FreqValidText;

			Wt::WContainerWidget *w;
			Wt::WBorderLayout *layout;

			std::vector<IonSpec> species;

			Wt::WGroupBox *West_Container, *ParameterType_Container, *IonSpecies_Container, 
					*PointParameterType_Container,*MagneticParameterType_Container;
			Wt::WButtonGroup *ParameterType_ButtonGroup;
			Wt::WComboBox *nIonSpecies_ComboBox;

			// Scan
			int nX;
			Wt::WText *nX_Text;
			Wt::WLineEdit *nX_LineEdit;

			// Magnetic field scan
			float r0, b0, rMin, rMax;
			Wt::WText *r0_Text, *b0_Text, *rMin_Text, *rMax_Text;
			Wt::WLineEdit *r0_LineEdit, *b0_LineEdit, *rMin_LineEdit, *rMax_LineEdit;

			Wt::Chart::WCartesianChart *b_CartesianChart;

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
			void AddIonSpecies();
			void RevealParameterSettings();

			int nSpecies;

};

HelloApplication::HelloApplication(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			w = new Wt::WContainerWidget(root());
			layout = new Wt::WBorderLayout();
			w->setLayout(layout,Wt::AlignTop | Wt::AlignJustify);

			//layout->addWidget(new Wt::WImage("naruto.jpg"),Wt::WBorderLayout::Center);

			enum CalculationType { Hot = 1, Cold = 2 };

			West_Container = new Wt::WGroupBox("Parameters");
			layout->addWidget(West_Container,Wt::WBorderLayout::West);
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
							boost::bind(&HelloApplication::RevealParameterSettings,this));

			PointParameterType_Container = new Wt::WGroupBox("Point Parameter Settings",root());

			Wt::WText *DensityText = new Wt::WText(PointParameterType_Container);
			DensityText->setText("Density [1/m^3]");
			DensityValidator = new Wt::WDoubleValidator(1e1,1e25,PointParameterType_Container);
			DensityValidator->setMandatory(1);
			Density_m3_LineEdit = new Wt::WLineEdit("2e19",PointParameterType_Container);
			Density_m3_LineEdit->setValidator(DensityValidator);
			DensityValidText = new Wt::WText(PointParameterType_Container);
			Density_m3_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,Density_m3_LineEdit,DensityValidText));

			new Wt::WBreak(PointParameterType_Container);

			Wt::WText *BFieldText = new Wt::WText(PointParameterType_Container);
			BFieldText->setText("BField [T]");
			BFieldValidator = new Wt::WDoubleValidator(1e-5,5e5,PointParameterType_Container);
			BFieldValidator->setMandatory(1);
			BField_LineEdit = new Wt::WLineEdit("5",PointParameterType_Container);
			BField_LineEdit->setValidator(BFieldValidator);
			BFieldValidText = new Wt::WText(PointParameterType_Container);
			BField_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,BField_LineEdit,BFieldValidText));

			new Wt::WBreak(PointParameterType_Container);

			Wt::WText *FreqText = new Wt::WText(PointParameterType_Container);
			FreqText->setText("Freq [Hz]");
			FreqValidator = new Wt::WDoubleValidator(1e1,5e12,PointParameterType_Container);
			FreqValidator->setMandatory(1);
			Freq_Hz_LineEdit = new Wt::WLineEdit("30e6",PointParameterType_Container);
			Freq_Hz_LineEdit->setValidator(FreqValidator);
			FreqValidText = new Wt::WText(PointParameterType_Container);
			Freq_Hz_LineEdit->changed().connect(
					boost::bind(&HelloApplication::editedLineEdit,this,Freq_Hz_LineEdit,FreqValidText));

			new Wt::WBreak(PointParameterType_Container);
			PointParameterType_Container->setHidden(1);

			MagneticParameterType_Container = new Wt::WGroupBox("Magnetic Field Scan Settings",root());
			MagneticParameterType_Container->setHidden(0);
			r0_Text = new Wt::WText("r0",MagneticParameterType_Container);
			r0_LineEdit = new Wt::WLineEdit("2",MagneticParameterType_Container);
			new Wt::WBreak(MagneticParameterType_Container);

			b0_Text = new Wt::WText("b0",MagneticParameterType_Container);
			b0_LineEdit = new Wt::WLineEdit("5",MagneticParameterType_Container);
			new Wt::WBreak(MagneticParameterType_Container);

			rMin_Text = new Wt::WText("rMin",MagneticParameterType_Container);
			rMin_LineEdit = new Wt::WLineEdit("0.1",MagneticParameterType_Container);
			new Wt::WBreak(MagneticParameterType_Container);

			rMax_Text = new Wt::WText("rMax",MagneticParameterType_Container);
			rMax_LineEdit = new Wt::WLineEdit("4",MagneticParameterType_Container);

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
				boost::bind(&HelloApplication::AddIonSpecies,this));

			new Wt::WBreak(IonSpecies_Container);

			nSpecies = 0;
			HelloApplication::AddIonSpecies();

			b_CartesianChart = new Wt::Chart::WCartesianChart();
			Wt::WStandardItemModel *model = new Wt::WStandardItemModel(40,2);
	  		model->setHeaderData(0, Wt::WString("X"));
	    	model->setHeaderData(1, Wt::WString("Y = sin(X)"));

		  	for (unsigned i = 0; i < 40; ++i) {

		      double x = (static_cast<double>(i) - 20) / 4;
		      model->setData(i, 0, x);
		      model->setData(i, 1, sin(x));

		    }

			b_CartesianChart->setModel(model);
			b_CartesianChart->setType(Wt::Chart::ScatterPlot);
			b_CartesianChart->setXSeriesColumn(0);
			b_CartesianChart->setLegendEnabled(true);
			b_CartesianChart->axis(Wt::Chart::XAxis).setLocation(Wt::Chart::ZeroValue);
			b_CartesianChart->resize(800,300);


			Wt::Chart::WDataSeries s(1,Wt::Chart::CurveSeries);
			b_CartesianChart->addSeries(s);

			layout->addWidget(b_CartesianChart,Wt::WBorderLayout::Center);


}

void HelloApplication::RevealParameterSettings()
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

void HelloApplication::AddIonSpecies()
{

		int nSpecTmp = nIonSpecies_ComboBox->currentIndex();
		std::cout << nSpecTmp << std::endl;

		for(int s=0; s<nSpecies; s++)
		{
				species[s].setHidden();
		}

		species.clear();
		nSpecies = 0;

		for(int s=0; s<=nSpecTmp; s++)
		{
				IonSpec mySpecParams(s+1,IonSpecies_Container);	
				species.push_back(mySpecParams);
				species[s].setVisible();
				nSpecies++;
		}

}

void HelloApplication::editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText)
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

void HelloApplication::greet()
{
		    greeting_->setText("Hello there, " + nameEdit_->text());
}

Wt::WApplication *createApplication(const Wt::WEnvironment& env)
{
		    return new HelloApplication(env);
}

int main(int argc, char **argv)
{
		    return Wt::WRun(argc, argv, &createApplication);
}
