import { Component, OnInit, Input, EventEmitter, Output } from '@angular/core';
import { FormBuilder, FormGroup, Validators } from '@angular/forms';
import { ParameterService, UserService } from '../../_service/index';
import { polarity, lc, statMethod, adjustP, species, instrument, ceInstt } from '../../_data_model/paras';

@Component({
  selector: 'app-analysis-para-setting',
  templateUrl: './analysis-para-setting.component.html',
  styleUrls: ['./analysis-para-setting.component.css']
})
export class AnalysisParaSettingComponent implements OnInit {
  public paraForm: FormGroup;
  // private paras: Paras;
  @Input() groups: string[];
  @Input() projectName: string;
  @Input() polarityMode: string;
  // private polarityMode: string = 'Positive';
  // private projectName: string = 'projectName2';
  // private groups: string[] = ['g1', 'g2'];
  private ceInstt = ceInstt;
  public polarityValue: string[] = [];
  public lcValue = lc;
  public insttValue = instrument;
  public ceValue = ceInstt['STT'];
  public statMethodValue = statMethod;
  public adjustPValue = adjustP;
  public speciesValue = species;
  private token: string = '';
  public visible = false;  // control Modal visible
  public visibleAnimate = false;  // control Modal visible
  public showError: boolean = false;
  public showSubmitLoader: boolean = false;
   @Output() finalSubmit = new EventEmitter<boolean>();

  constructor(
    private fb: FormBuilder,
    private paraService: ParameterService,
    private userService: UserService
  ) { }

  ngOnInit() {
    this.polarityValue = polarity[this.polarityMode];
    let currentUser = this.userService.getCurrentUser();
    if (currentUser) {
      this.token = currentUser['token'];
    }
    // console.log(this.groups);
    // console.log(this.projectName);
    if (this.polarityValue.length !== 0) {
      this.createForm();
    }

    this.paraForm.controls['instrument'].valueChanges.subscribe((value) => {
      // console.log(value);
      if (value === 'Sciex TripleTOF') {
        this.ceValue = this.ceInstt['STT'];
      }
      else if (value === 'Agilent QTOF') {
        this.ceValue = this.ceInstt['AQT'];
      }
      else if (value === 'Other QTOF') {
        this.ceValue = this.ceInstt['OQ'];
      }
      else if (value === 'Thermo Orbitrap (HCD)') {
        this.ceValue = this.ceInstt['TOQ'];
      }
      this.paraForm.controls['ce'].setValue(this.ceValue[0]);
    })
  }

  createForm() {
    this.paraForm = this.fb.group({
      polarity: [this.polarityValue[0], Validators.required],
      lcType: [this.lcValue[0], Validators.required],
      instrument: [this.insttValue[0], Validators.required],
      ce: [this.ceValue[0], Validators.required],
      controlGroup: [this.groups[0], Validators.required],
      caseGroup: [this.groups[1], Validators.required],
      statMethod: [this.statMethodValue[0], Validators.required],
      adjustP: [this.adjustPValue[0], Validators.required],
      pCutoff: [0.05, Validators.required],
      species: [this.speciesValue[0], Validators.required]
    });
  }

  onSubmit() {
    // console.log(this.paraForm);
    // console.log(typeof(this.paraForm.value));
    // console.log(this.token);
    // this.show();
    this.showSubmitLoader = true;
    
    let createQueueInfo = this.paraService.createProjectQueue(this.token, this.projectName, 
                                                              this.paraForm.value);
    setTimeout(() => {
      createQueueInfo
      .then(res => {
        // console.log(res);
        if (res===201) {
          this.show();
          this.showSubmitLoader = false;
        }
      })
      .catch(res => {
        this.showError = true;
        this.showSubmitLoader = false;
        // console.log(res);
    });
    }, 1);

  }

  public show(): void {
    this.visible = true;
    setTimeout(() => this.visibleAnimate = true, 100);
  }

  public hide(): void {
    this.visibleAnimate = false;
    setTimeout(() => this.visible = false, 100);
    this.finalSubmit.emit(true);
    this.paraForm.reset();
  }

}
