import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { AnalysisParaSettingComponent } from './analysis-para-setting.component';

describe('AnalysisParaSettingComponent', () => {
  let component: AnalysisParaSettingComponent;
  let fixture: ComponentFixture<AnalysisParaSettingComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ AnalysisParaSettingComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(AnalysisParaSettingComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should be created', () => {
    expect(component).toBeTruthy();
  });
});
