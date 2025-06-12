import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Open Source',
    Svg: require('@site/static/img/open_source.svg').default,
    description: (
      <>
        TidyScreen was originally developed and is permanently mantained with focus 
	on the integration of Open Source tools.
      </>
    ),
  },
  {
    title: 'Reproducible screening workflows',
    Svg: require('@site/static/img/reproducibility.svg').default,
    description: (
      <>
        The main objective of TidyScreen is to turn complex <i>in silico</i> screening pipelines into highly reproducible workflows.
      </>
    ),
  },
  {
    title: 'Community driven development',
    Svg: require('@site/static/img/feedback.svg').default,
    description: (
      <>
        Each <i>in silico</i> drug screening campaing encompasses its own particularities, and TidyScreen is expected to grow upon use cases feedback.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} role="img" />
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
