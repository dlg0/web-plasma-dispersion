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

class IonSpecParameters
{
		public:
		IonSpecParameters(void);
		~IonSpecParameters(void);

		private:

		public:
		Wt::WLineEdit *SpecAMU_LineEdit, *SpecZ_LineEdit;
		Wt::WBreak *SpecBreak;
		Wt::WText *SpecText;
};

IonSpecParameters::IonSpecParameters(void)
{
}
IonSpecParameters::~IonSpecParameters(void)
{
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

			std::vector<IonSpecParameters> species;

			Wt::WLineEdit *Spec1AMU, *Spec2AMU, *Spec3AMU, *Spec4AMU, *Spec5AMU, *Spec6AMU;
			Wt::WLineEdit *Spec1Z, *Spec2Z, *Spec3Z, *Spec4Z, *Spec5Z, *Spec6Z;
			Wt::WText *Spec1Text, *Spec2Text, *Spec3Text, *Spec4Text, *Spec5Text, *Spec6Text;
			Wt::WBreak *Spec1Break, *Spec2Break, *Spec3Break, *Spec4Break, *Spec5Break, *Spec6Break;
			Wt::WGroupBox *HotOrCold_Container, *ParameterType_Container, *IonSpecies_Container, 
					*PointParameterType_Container,*MagneticParameterType_Container;
			Wt::WButtonGroup *ParameterType_ButtonGroup;
			Wt::WComboBox *nIonSpecies_ComboBox;

			void greet();
			void editedLineEdit(Wt::WLineEdit *LineEditIn, Wt::WText *ValidText);
			void showIonSpecies();
			void RevealParameterSettings();

};

HelloApplication::HelloApplication(const Wt::WEnvironment& env) : Wt::WApplication(env)
{
		    setTitle("Plasma Dispersion Calculator");

			//root()->addWidget(new Wt::WText("Your name, please ? "));
			//nameEdit_ = new Wt::WLineEdit(root());
			//Wt::WPushButton *button = new Wt::WPushButton("RUN", root());
			//root()->addWidget(new Wt::WBreak());
			//greeting_ = new Wt::WText(root());
			//button->clicked().connect(this, &HelloApplication::greet);

			enum CalculationType { Hot = 1, Cold = 2 };

			HotOrCold_Container = new Wt::WGroupBox("Hot or Cold Calculation",root());
			Wt::WButtonGroup *HotOrCold_ButtonGroup = new Wt::WButtonGroup(root());

			Wt::WRadioButton *HotOrCold_Button;

			HotOrCold_Button = new Wt::WRadioButton("Hot", HotOrCold_Container);
			new Wt::WBreak(HotOrCold_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Hot);

			HotOrCold_Button = new Wt::WRadioButton("Cold", HotOrCold_Container);
			new Wt::WBreak(HotOrCold_Container);
			HotOrCold_ButtonGroup->addButton(HotOrCold_Button,Cold);

			HotOrCold_ButtonGroup->setCheckedButton(HotOrCold_ButtonGroup->button(Cold));

			enum ParameterType { Point = 1, MagneticField = 2, Numerical = 3 };

			ParameterType_Container = new Wt::WGroupBox("Parameter Scan Type",root());
			ParameterType_ButtonGroup = new Wt::WButtonGroup(root());

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
				boost::bind(&HelloApplication::showIonSpecies,this));

			new Wt::WBreak(IonSpecies_Container);

			Spec1Text = new Wt::WText("Ion Species 1 (AMU, Z)",IonSpecies_Container);
			Spec1AMU = new Wt::WLineEdit("2",IonSpecies_Container);
			Spec1Z = new Wt::WLineEdit("1",IonSpecies_Container);
			Spec1Break = new Wt::WBreak(IonSpecies_Container);

			Spec2Text = new Wt::WText("Ion Species 2 (AMU, Z)",IonSpecies_Container);
			Spec2Text->setHidden(1);
			Spec2AMU = new Wt::WLineEdit(IonSpecies_Container);
			Spec2AMU->setHidden(1);
			Spec2Z = new Wt::WLineEdit(IonSpecies_Container);
			Spec2Z->setHidden(1);
			Spec2Break = new Wt::WBreak(IonSpecies_Container);
			Spec2Break->setHidden(1);

			Spec3Text = new Wt::WText("Ion Species 3 (AMU, Z)",IonSpecies_Container);
			Spec3Text->setHidden(1);
			Spec3AMU = new Wt::WLineEdit(IonSpecies_Container);
			Spec3AMU->setHidden(1);
			Spec3Z = new Wt::WLineEdit(IonSpecies_Container);
			Spec3Z->setHidden(1);
			Spec3Break = new Wt::WBreak(IonSpecies_Container);
			Spec3Break->setHidden(1);

			Spec4Text = new Wt::WText("Ion Species 4 (AMU, Z)",IonSpecies_Container);
			Spec4Text->setHidden(1);
			Spec4AMU = new Wt::WLineEdit(IonSpecies_Container);
			Spec4AMU->setHidden(1);
			Spec4Z = new Wt::WLineEdit(IonSpecies_Container);
			Spec4Z->setHidden(1);
			Spec4Break = new Wt::WBreak(IonSpecies_Container);
			Spec4Break->setHidden(1);

			Spec5Text = new Wt::WText("Ion Species 5 (AMU, Z)",IonSpecies_Container);
			Spec5Text->setHidden(1);
			Spec5AMU = new Wt::WLineEdit(IonSpecies_Container);
			Spec5AMU->setHidden(1);
			Spec5Z = new Wt::WLineEdit(IonSpecies_Container);
			Spec5Z->setHidden(1);
			Spec5Break = new Wt::WBreak(IonSpecies_Container);
			Spec5Break->setHidden(1);

			Spec6Text = new Wt::WText("Ion Species 6 (AMU, Z)",IonSpecies_Container);
			Spec6Text->setHidden(1);
			Spec6AMU = new Wt::WLineEdit(IonSpecies_Container);
			Spec6AMU->setHidden(1);
			Spec6Z = new Wt::WLineEdit(IonSpecies_Container);
			Spec6Z->setHidden(1);
			Spec6Break = new Wt::WBreak(IonSpecies_Container);
			Spec6Break->setHidden(1);

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

void HelloApplication::showIonSpecies()
{

		int nSpec = nIonSpecies_ComboBox->currentIndex();
		std::cout << nSpec << std::endl;
		nSpec++;

		IonSpecParameters mySpecParams;	
		for(int s=0; s<nSpec; s++)
		{
				//Fill in mySpecParams then add to the species vec.
				// ....
				species.push_back(mySpecParams);
		}

		if(nSpec>=1){
			Spec1Text->setHidden(0);
			Spec1AMU->setHidden(0);
			Spec1Z->setHidden(0);
			Spec1Break->setHidden(0);
		}
		else
		{
			Spec1Text->setHidden(1);
			Spec1AMU->setHidden(1);
			Spec1Z->setHidden(1);
			Spec1Break->setHidden(1);
		}

		if(nSpec>=2){
			Spec2Text->setHidden(0);
			Spec2AMU->setHidden(0);
			Spec2Z->setHidden(0);
			Spec2Break->setHidden(0);
		}
		else
		{
			Spec2Text->setHidden(1);
			Spec2AMU->setHidden(1);
			Spec2Z->setHidden(1);
			Spec2Break->setHidden(1);
		}
	
		if(nSpec>=3){
			Spec3Text->setHidden(0);
			Spec3AMU->setHidden(0);
			Spec3Z->setHidden(0);
			Spec3Break->setHidden(0);
		}
		else
		{
			Spec3Text->setHidden(1);
			Spec3AMU->setHidden(1);
			Spec3Z->setHidden(1);
			Spec3Break->setHidden(1);
		}
	
		if(nSpec>=4){
			Spec4Text->setHidden(0);
			Spec4AMU->setHidden(0);
			Spec4Z->setHidden(0);
			Spec4Break->setHidden(0);
		}
		else
		{
			Spec4Text->setHidden(1);
			Spec4AMU->setHidden(1);
			Spec4Z->setHidden(1);
			Spec4Break->setHidden(1);
		}
		
		if(nSpec>=5){
			Spec5Text->setHidden(0);
			Spec5AMU->setHidden(0);
			Spec5Z->setHidden(0);
			Spec5Break->setHidden(0);
		}
		else
		{
			Spec5Text->setHidden(1);
			Spec5AMU->setHidden(1);
			Spec5Z->setHidden(1);
			Spec5Break->setHidden(1);
		}

		if(nSpec>=6){
			Spec6Text->setHidden(0);
			Spec6AMU->setHidden(0);
			Spec6Z->setHidden(0);
			Spec6Break->setHidden(0);
		}
		else
		{
			Spec6Text->setHidden(1);
			Spec6AMU->setHidden(1);
			Spec6Z->setHidden(1);
			Spec6Break->setHidden(1);
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
