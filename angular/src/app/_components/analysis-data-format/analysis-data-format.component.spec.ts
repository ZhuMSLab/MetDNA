import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { AnalysisDataFormatComponent } from './analysis-data-format.component';

describe('AnalysisDataFormatComponent', () => {
  let component: AnalysisDataFormatComponent;
  let fixture: ComponentFixture<AnalysisDataFormatComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ AnalysisDataFormatComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(AnalysisDataFormatComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should be created', () => {
    expect(component).toBeTruthy();
  });
});
